import argparse
import re
import shutil
import subprocess
import sys
from typing import List, Optional, Tuple


def _run_command(cmd):
    return subprocess.run(cmd, check=False, capture_output=True, text=True)


def _find_conda_executable() -> Optional[str]:
    return shutil.which("conda")


def _detect_cuda_version() -> Tuple[Optional[int], Optional[int], str]:
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return None, None, "nvidia-smi not found"

    result = _run_command([nvidia_smi])
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or "nvidia-smi failed"
        return None, None, message

    match = re.search(r"CUDA Version:\s*([0-9]+)\.([0-9]+)", result.stdout)
    if not match:
        return None, None, "CUDA version not found in nvidia-smi output"

    return int(match.group(1)), int(match.group(2)), "ok"


def _select_cudatoolkit_pin(cuda_major: int, cuda_minor: int) -> Tuple[Optional[str], str]:
    if cuda_major < 11:
        return None, "Detected driver CUDA capability is below 11.x; skipping CUDA toolkit pin and recommending CPU platform."

    if cuda_major == 11:
        return f"cudatoolkit={cuda_major}.{cuda_minor}", f"Detected CUDA {cuda_major}.{cuda_minor}; pinning matching cudatoolkit."

    return "cudatoolkit=11.8", (
        f"Detected CUDA {cuda_major}.{cuda_minor}. "
        "Using cudatoolkit=11.8 for broad OpenMM/conda compatibility on modern drivers."
    )


def _build_install_command(conda_exe: str, env_name: str, channels, package_specs):
    cmd = [
        conda_exe,
        "install",
        "-n",
        env_name,
        "--override-channels",
        "-y",
    ]
    for channel in channels:
        cmd.extend(["-c", channel])
    cmd.extend(package_specs)
    return cmd


def _run_openmm_installation_test(conda_exe: str, env_name: str):
    return _run_command([conda_exe, "run", "-n", env_name, "python", "-m", "openmm.testInstallation"])


def _is_ptx_or_cuda_runtime_error(test_result) -> bool:
    combined = f"{test_result.stdout}\n{test_result.stderr}"
    lowered = combined.lower()
    return (
        "cuda_error_unsupported_ptx_version" in lowered
        or "error loading cuda module" in lowered
        or "unsupported ptx" in lowered
    )


def _build_retry_pins(cuda_major: Optional[int], cuda_minor: Optional[int], initial_pin: Optional[str]) -> List[str]:
    pins = []
    if initial_pin:
        pins.append(initial_pin)

    if cuda_major == 11 and cuda_minor is not None:
        for mn in range(cuda_minor, 3, -1):
            pins.append(f"cudatoolkit=11.{mn}")
    else:
        pins.extend(["cudatoolkit=11.8", "cudatoolkit=11.7", "cudatoolkit=11.6", "cudatoolkit=11.4"])

    deduped = []
    for pin in pins:
        if pin and pin not in deduped:
            deduped.append(pin)
    return deduped


def main() -> int:
    parser = argparse.ArgumentParser(description="Cross-platform bootstrap installer for mdpertool conda environment.")
    parser.add_argument("--env-name", default="mdpertool", help="Conda environment name")
    parser.add_argument("--python-version", default="3.9", help="Python version for the environment")
    parser.add_argument("--cpu-only", action="store_true", help="Skip CUDA toolkit pinning and install CPU-oriented setup")
    parser.add_argument("--cudatoolkit-pin", default=None, help="Override automatic toolkit pin, e.g. cudatoolkit=11.8")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument(
        "--channels",
        nargs="+",
        default=["bio-otto", "conda-forge"],
        help="Conda channel priority order",
    )

    args = parser.parse_args()

    conda_exe = _find_conda_executable()
    if not conda_exe:
        print("ERROR: conda executable was not found in PATH.")
        print("Please open an Anaconda/Miniconda shell or initialize conda in your shell.")
        return 1

    create_cmd = [
        conda_exe,
        "create",
        "-n",
        args.env_name,
        f"python={args.python_version}",
        "-y",
    ]

    package_specs = ["mdpertool"]
    detected_major = None
    detected_minor = None
    auto_selected_pin = None

    if args.cpu_only:
        print("INFO: CPU-only mode enabled; no CUDA toolkit pin will be applied.")
    else:
        if args.cudatoolkit_pin:
            package_specs.append(args.cudatoolkit_pin)
            print(f"INFO: Using user-provided toolkit pin: {args.cudatoolkit_pin}")
        else:
            major, minor, reason = _detect_cuda_version()
            detected_major, detected_minor = major, minor
            if major is None or minor is None:
                print(f"INFO: GPU/CUDA auto-detection unavailable ({reason}). Installing without toolkit pin.")
            else:
                pin, pin_reason = _select_cudatoolkit_pin(major, minor)
                print(f"INFO: {pin_reason}")
                if pin:
                    auto_selected_pin = pin
                    package_specs.append(pin)

    install_cmd = _build_install_command(conda_exe, args.env_name, args.channels, package_specs)

    print("\nPlanned commands:")
    print(" ".join(create_cmd))
    print(" ".join(install_cmd))

    if args.dry_run:
        print("\nDry-run mode: no commands were executed.")
        return 0

    create_result = _run_command(create_cmd)
    if create_result.returncode != 0:
        print("\nERROR: Failed to create conda environment.")
        print(create_result.stdout)
        print(create_result.stderr)
        return create_result.returncode

    install_result = _run_command(install_cmd)
    if install_result.returncode != 0:
        print("\nERROR: Failed to install mdpertool packages.")
        print(install_result.stdout)
        print(install_result.stderr)
        return install_result.returncode

    if not args.cpu_only:
        test_result = _run_openmm_installation_test(conda_exe, args.env_name)
        if test_result.returncode == 0:
            print("\nINFO: OpenMM installation test passed.")
        elif _is_ptx_or_cuda_runtime_error(test_result):
            print("\nWARNING: CUDA/PTX compatibility issue detected. Trying automatic toolkit pin fallback...")
            retry_pins = _build_retry_pins(detected_major, detected_minor, args.cudatoolkit_pin or auto_selected_pin)

            recovered = False
            for pin in retry_pins:
                print(f"INFO: Retrying with {pin}")
                retry_cmd = _build_install_command(conda_exe, args.env_name, args.channels, ["mdpertool", pin])
                retry_install = _run_command(retry_cmd)
                if retry_install.returncode != 0:
                    continue

                retry_test = _run_openmm_installation_test(conda_exe, args.env_name)
                if retry_test.returncode == 0:
                    recovered = True
                    print(f"INFO: CUDA/OpenMM test passed after pinning {pin}")
                    break

            if not recovered:
                print("WARNING: Automatic CUDA toolkit fallback did not fully resolve CUDA runtime mismatch.")
                print("WARNING: MDPerTool environment is installed; set platform to CPU in the app until driver/toolkit are aligned.")
        else:
            print("WARNING: OpenMM installation test reported issues not recognized as CUDA PTX mismatch.")
            print(test_result.stdout)
            print(test_result.stderr)

    print("\nSUCCESS: mdpertool environment is ready.")
    print(f"Activate with: conda activate {args.env_name}")
    print("Run app: mdpertool gui")
    return 0


if __name__ == "__main__":
    sys.exit(main())

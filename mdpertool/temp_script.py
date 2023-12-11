import re
msg = "INFO: Res Num: 11, Res Num: 447, Decomposition Progress: 2.470850111856823"


if 'Decomposition Progress:' in msg:
    decomp_info_log = re.search(r"Decomposition Progress: ([\d.]+)", msg)

    if decomp_info_log:
        decompose_started = True
        extracted_number = float(decomp_info_log.group(1))
        formatted_number = round(extracted_number, 2)
        print(formatted_number)
import openmm as mm
import openmm.app as app

def test_platform(name, precision='single'):
    print(f"\nTesting {name} with precision: {precision}")
    try:
        platform = mm.Platform.getPlatformByName(name)
        # Create a dummy system to test context creation
        system = mm.System()
        integrator = mm.VerletIntegrator(0.001)
        properties = {}
        if name in ['OpenCL', 'CUDA']:
            properties[f"{name}Precision"] = precision
        
        # We need a topology and positions for a real simulation, 
        # but just creating a Context might be enough to trigger the error.
        # Actually, let's just try to create a context if possible.
        # context = mm.Context(system, integrator, platform, properties)
        # Wait, Context needs at least one particle.
        system.addParticle(1.0)
        context = mm.Context(system, integrator, platform, properties)
        print(f"  SUCCESS: Context created with {name} / {precision}")
        del context
    except Exception as e:
        print(f"  FAILURE: {e}")

test_platform('OpenCL', 'single')
test_platform('OpenCL', 'mixed')
test_platform('OpenCL', 'double')

test_platform('CPU')
test_platform('Reference')

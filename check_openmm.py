import openmm as mm
print("Available Platforms:")
for i in range(mm.Platform.getNumPlatforms()):
    p = mm.Platform.getPlatform(i)
    print(f"Platform {i}: {p.getName()}")
    try:
        if p.getName() in ['OpenCL', 'CUDA']:
             print(f"  Precision levels: {p.getPropertyValue(None, 'Precision') if 'Precision' in p.getPropertyNames() else 'Unknown'}")
    except:
        pass

print("\nDefault Platform:", mm.Platform.getPlatform(0).getName())

try:
    platform = mm.Platform.getPlatformByName('OpenCL')
    print("\nOpenCL Platform found.")
    print("OpenCL Properties:", platform.getPropertyNames())
except Exception as e:
    print("\nOpenCL Platform NOT found:", e)

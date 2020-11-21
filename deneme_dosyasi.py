import os

def parseFloat(str):
    try:
        return float(str)
    except:
        str = str.strip()
        if str.endswith("%"):
            return float(str.strip("%").strip()) / 100
        raise Exception("Don't know how to parse %s" % str)




text = """#"Progress (%)","Step","Potential Energy (kJ/mole)","Temperature (K)","Speed (ns/day)","Time Remaining"
100.0%,1000,-3478383.1367646214,277.9379980064646,--,--"""

initial, line = text.split('\n', 1)
headers = initial.strip().split(',')
_headers = [e.strip('#"\'') for e in headers]
t =[e.strip('%"\'') for e in line.strip().split(',')]
print(t[1])
# values = map(float, t)
msg = dict(zip(_headers, t))
print(msg)

# headers = initial.strip().split(self._separator)
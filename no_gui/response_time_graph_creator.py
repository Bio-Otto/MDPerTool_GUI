import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from ast import literal_eval

infile = r"api.log"

keep_phrases = [
    "Topology File:",
    "Speed Factor:",
    "Perturbed Residue(s):"
]

top_file_name = None
speed_factors = None
perturbed_residues = None

with open(infile) as f:
    f = f.readlines()

for line in f:
    for phrase in keep_phrases:
        if phrase in line:
            if phrase == "Topology File:":
                top_file_name = os.path.basename(line.strip().split('\t')[1]).split('.')[0]

            if phrase == "Speed Factor:":
                speed_factors = np.array(literal_eval(line.strip().split('\t')[1]))

            if phrase == "Perturbed Residue(s):":
                perturbed_residues = '-'.join(np.array(literal_eval(line.strip().split('\t')[1])))

            break


def getResponseTimeGraph(responseTimeFile, factor):
    # Draw Graph residue response vs Time
    Responses_file = pd.read_csv(responseTimeFile, header=None)
    Response_Time_Column = Responses_file.values
    Response_Time_Column_Max = Response_Time_Column.max()
    row, col = Responses_file.shape
    Response_Count = []
    Increase = 0

    for i in range(int(Response_Time_Column_Max)):
        Increase = Increase + list(Response_Time_Column).count(i)
        Response_Count.append(Increase)
    plt.ylim(0, row + 5)
    plt.plot(Response_Count, label='x%s' % factor, linewidth=1)

    return row, col


if speed_factors is not None:
    for factor in sorted(speed_factors):
        try:
            row_num, col_num = getResponseTimeGraph('responseTimes_%s.csv' % factor, factor=factor)
        except:
            print("There is no % named file in the directory")

plt.xlabel('Time Step (Frame)')
plt.ylabel('Responded Residue Count')
plt.title('Response Time Graph: %s' % top_file_name)

ax = plt.gca() # get axis handle
max_handler = 0
for i, line in enumerate(ax.lines):
    x_data = ax.lines[i].get_xdata() # get the first line, there might be more
    try:
        max_val = max(x_data)
        if max_val > max_handler:
            max_handler = max_val

    except Exception as ERr:
        print(ERr)

plt.text(max_handler-350, 10, 'Perturbed Residue(s): %s' % perturbed_residues, style='italic',
         bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})
plt.legend(title='Speed Factor')

plt.margins(x=0.01, y=0.01, tight=True)
plt.tight_layout()

fig = plt.gcf()
plt.show()  # show it here (important, if done before you will get blank picture)
#fig.set_size_inches((25, 15), forward=False)
fig.savefig('%s.png' % top_file_name, dpi=300)  # Change is over here

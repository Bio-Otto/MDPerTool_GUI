import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def getResponseTimeGraph(responseTimeFile):
    # Draw Graph residue response vs Time

    Responses_file = pd.read_csv(responseTimeFile, header=None)
    Response_Time_Column = Responses_file.values
    Response_Time_Column_Max = Response_Time_Column.max()

    Response_Count = []
    Increase = 0

    for i in range(int(Response_Time_Column_Max)):
        Increase = Increase + list(Response_Time_Column).count(i)
        Response_Count.append(Increase)

    plt.plot(Response_Count)



rfile = 'responseTimes.csv'
getResponseTimeGraph(rfile)
# r2 = 'responseTimes_amber2000.csv'
# getResponseTimeGraph(r2)
plt.show()
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def getCoordsForPlot(step):
    df = pd.read_csv("24thousand.csv")
    df = df[["Step", "Total cells", "Position", "Cell Type"]]

    # Select the step you want to plot, from 0 to 24000 (~11 days)
    df_step0 = df.loc[df["Step"] == step]

    mPoints = df_step0.loc[df_step0["Cell Type"] == "mesenchymal"]["Position"] # Series object
    ePoints = df_step0.loc[df_step0["Cell Type"] == "epithelial"]["Position"]
    #vPoints = df_step0.loc[(df_step0["Agent Type"] == "vessel") & (df_step0["Ruptured"] == False)]
    #vRupturedPoints = df_step0.loc[(df_step0["Agent Type"] == "vessel") & (df_step0["Ruptured"] == True)]

    mPoints = list(mPoints.map(eval)) # [(104, 101), (101, 97), (101, 95)]
    ePoints = list(ePoints.map(eval))

    Xm, Ym = [i[0] for i in mPoints], [i[1] for i in mPoints]
    Xe, Ye = [i[0] for i in ePoints], [i[1] for i in ePoints]

    return [Xm, Ym, Xe, Ye]

def plotGrowthData(fig_index):
    # get data at each step
    df = pd.read_csv("24thousand.csv")
    df = df[["Step", "Total cells", "Cell Type"]]

    plt.figure(fig_index)

    # For mesenchymal
    stepsize = 3000 #doubling rate
    df_m = df.loc[(df["Step"] % stepsize == 0)]
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_m["Step"]) + 1, stepsize)
    numberMesenchymalEachStep = [df_m.loc[(df_m["Step"] == step) & (df_m["Cell Type"] == "mesenchymal")].shape[0] for step in steps]
    plt.plot(steps*11/24000, numberMesenchymalEachStep, label="Mesenchymal cells")

    # For epithelial
    stepsize = 2000 #doubling rate
    df_e = df.loc[(df["Step"] % stepsize == 0)]
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_e["Step"]) + 1, stepsize)
    numberEpithelialEachStep = [df_e.loc[(df_e["Step"] == step) & (df_e["Cell Type"] == "epithelial")].shape[0] for step in steps]
    
    plt.plot(steps*11/24000, numberEpithelialEachStep, label="Epithelial cells")

    plt.title('Tumor growth')
    plt.xlabel('Days')
    plt.ylabel('Number of cells')
    plt.legend(loc="upper left")

    plt.show()


def plotCancer(coordsList,i):

    Xm, Ym, Xe, Ye = coordsList[0], coordsList[1], coordsList[2], coordsList[3]
    plt.figure(i, figsize=(6, 5))
    
    plt.scatter(Xm, Ym, marker='o', alpha=0.10, label="Mesenchymal cells")
    plt.scatter(Xe, Ye, marker='h', alpha=0.05, label="Epithelial cells")
    plt.xlim(0, 201)
    plt.ylim(0, 201)

    xticks = np.arange(0, 201, step=40) # 6 ticks
    xticklabels = [str(round(i,1)) for i in np.arange(0, 2.1, step = 2/201*40)]
    plt.xticks(xticks, xticklabels)
    plt.yticks(xticks, xticklabels)
    plt.xlabel("mm")
    plt.ylabel("mm")
    
    #plt.style.use("seaborn")

if __name__ == "__main__":
    
    # Plot a certain state of the metastasis
    # t0 : step = 0
    # 11 days: step = 24000 
    step = 0 

    for id, step in enumerate(reversed(range(0,24000+1,2000))):
        plotCancer(getCoordsForPlot(step),id)
        if step == 0:
            plt.title(f'Initial Tumor')
            plt.legend(loc="upper left")
        else:
            plt.title(f'Tumor size at {11/24000 * step:.2f} days ({step} steps)')
    
    plotGrowthData(id+1)
      
    
    plt.show()






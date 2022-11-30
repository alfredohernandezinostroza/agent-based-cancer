import pandas as pd
import matplotlib.pyplot as plt

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


def plotCancer(coordsList,i):

    Xm, Ym, Xe, Ye = coordsList[0], coordsList[1], coordsList[2], coordsList[3]
    plt.figure(i, figsize=(6, 5))
    
    plt.scatter(Xm, Ym, marker='o', alpha=0.25)
    plt.scatter(Xe, Ye, marker='.', alpha=0.15)
    plt.xlim(0, 201)
    plt.ylim(0, 201)
 

if __name__ == "__main__":
    
    # Plot a certain state of the metastasis
    # t0 : step = 0
    # 11 days: step = 24000 
    step = 0 

    for id, step in enumerate(reversed(range(0,24000+1,2000))):
        plotCancer(getCoordsForPlot(step),id)
        plt.title(f'Tumor metastasys at {11/24000 * step:.2f} days (step {step}')
    
    plt.show()






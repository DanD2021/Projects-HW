import pandas as pd
import numpy as np
import sys
import mykmeanssp as kmean

if len(sys.argv) != 4:
    print("Invalid Input!")
    sys.exit()

k = int(sys.argv[1])

if k == 1:
    print("Invalid Input!")
    sys.exit()

if sys.argv[2] == "spk":
    goal = 0

elif sys.argv[2] == "wam":
    goal = 1

elif sys.argv[2] == "ddg":
    goal = 2

elif sys.argv[2] == "lnorm":
    goal = 3

elif sys.argv[2] == "jacobi":
    goal = 4

else:
    print("Invalid Input!")
    sys.exit()

df = pd.read_csv(sys.argv[3], header=None)
numberOfPoints, pointSize = df.shape
listToC = df.to_numpy().flatten()
listToC = list(listToC[~np.isnan(listToC)])

if (sys.argv[2] == "spk") and (numberOfPoints < k):
    print("Invalid Input!")
    sys.exit()
    
list_from_c = kmean.new_fit(listToC, numberOfPoints, pointSize, k, goal)

if list_from_c == None:
    print("An Error Has Occurred")
    sys.exit()
if list_from_c == NameError:
    print("An Error Has Occurred")
    sys.exit()
if list_from_c == 0:
    print("An Error Has Occurred")
    sys.exit()

pointSize = int(list_from_c[-2])

if 0 < goal:
    
    if goal==4:
        numberOfPoints+=1
    
    for i in range(numberOfPoints):
        clus = ""
        for j in range(pointSize):
            clus += format(list_from_c[i * pointSize + j], '.4f')
            if j != pointSize - 1:
                clus += ","
        print(clus)

else:
    if k == 0:
        k = int(list_from_c[-1])
        if k==1:
            print("An Error Has Occurred")
            sys.exit()
    df = []

    for i in range(numberOfPoints):
        df.append([])
        for j in range(pointSize):
            df[i].append(list_from_c[i*pointSize + j])

    df = pd.DataFrame(df)
    matOfPoints = df.to_numpy()  # convert df to numpy instance

    np.random.seed(0)
    randomPoint = np.random.choice(numberOfPoints)  # choose random index

    clusters = []
    test = []

    test.append(randomPoint)    # add to test  the random choosen index
    clusters.append(matOfPoints[randomPoint])   # add the point to the new numpy  with the random index

    distance = [-1 for i in range(numberOfPoints)]
    pro = [0 for i in range(numberOfPoints)]

    while len(clusters) < k:
        i = 0
        for point in matOfPoints:
            distance[i] = -1
            for cluster in clusters:
                dotPoint = point - cluster
                dotP = dotPoint.dot(dotPoint)
                if distance[i] == -1:
                    distance[i] = dotP
                else:
                    distance[i] = min(dotP, distance[i])

            i += 1

        sumD = sum(distance)
        for i in range(numberOfPoints):
            pro[i] = distance[i] / sumD
        x = np.random.choice(numberOfPoints, p=pro)
        test.append(x)
        clusters.append(matOfPoints[x])

    print(','.join(str(round(x, 4)) for x in test))

    listToC = matOfPoints.tolist()
    listOfPoints = [i.tolist() for i in clusters]
    for index in range(len(listToC)):
        if index not in test:
            listOfPoints.append(listToC[index])
    listToC = [i for coord in listOfPoints for i in coord]

    numberOfClusters = k
    lenOflist = (k * pointSize)

    listToPrint = kmean.fit(listToC, numberOfPoints, pointSize, numberOfClusters, 300, 0.0, lenOflist)
    
    if listToPrint is None:
        sys.exit()
    if listToPrint == NameError:
        sys.exit()
    if listToPrint == 0:
        sys.exit()

    for i in range(numberOfClusters):
        clus = ""
        for j in range(pointSize):
            clus += format(listToPrint[i * pointSize + j], '.4f')
            if j != pointSize - 1:
                clus += ","
        print(clus)

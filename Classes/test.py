import random

def disaggregate_clusters(vasculature, time):
    big_clusters = [cluster for cluster in vasculature[time] if sum(cluster) > 1]
    new_vasculature = [cluster for cluster in vasculature[time] if sum(cluster) == 1]
    for cluster in big_clusters:
        new_mesenchymal, new_epithelial = cluster
        for ccell_type, ccells_amount in enumerate(cluster):
            for i in range(ccells_amount):
                if random.random() > 0.5:
                    if ccell_type == 0:
                        new_vasculature += [(1, 0)]
                        new_mesenchymal -= 1
                    if ccell_type == 1:
                        new_vasculature += [(0, 1)]
                        new_epithelial -= 1
        if new_mesenchymal + new_epithelial > 0:
            new_vasculature += [(new_mesenchymal,new_epithelial)]
    vasculature[time] = new_vasculature
    print(new_vasculature)


                
vasculature = { 4: [(1,1)],
                5: [(3,4), (2,2), (1,1)]}

test = 5
print(sum([sum(t) for t in vasculature[test]]))
disaggregate_clusters(vasculature, test)
print(sum([sum(t) for t in vasculature[test]]))

import random


class Order():

    def __init__(self, km1Size):
        self.km1Size = km1Size



    def forwardSort(self, query_str, mapF, ID, i, fsD):

        lastx = query_str[-self.km1Size:]

        map = mapF[ID]

        if lastx in map.keys():
            possAdds = map[lastx]

        else:
            possAdds = []

        if len(possAdds) > 0:
            add = random.choice(possAdds)
            map[lastx].remove(add)
            newQuery = "".join([query_str, add[-1:]])

            if ID not in fsD.keys():
                fsD[ID] = [i]

            else:
                new = []
                for n in fsD[ID]:
                    new.append(n)

                new.append(i)
                fsD[ID] = new

            return newQuery

        else:

            return query_str

    def backwardSort(self, query_str, mapB, ID, i, bsD):
        firstx = query_str[:self.km1Size]

        map = mapB[ID]

        if firstx in map.keys():
            possAdds = map[firstx]

        else:
            possAdds = []

        if len(possAdds) > 0:
            add = random.choice(possAdds)
            map[firstx].remove(add)
            newQuery = "".join([add[:1], query_str])

            if ID not in bsD.keys():
                bsD[ID] = [i]
            else:
                new = []
                for n in bsD[ID]:
                    new.append(n)

                new.append(i)
                bsD[ID] = new

            return newQuery

        else:

            return query_str


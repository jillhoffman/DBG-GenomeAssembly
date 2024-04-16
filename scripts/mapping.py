class map():
    def __init__(self, map):
        self.km1Map = map


    def makeMaps(self, start, target):

        if start not in self.km1Map.keys():
            self.km1Map[start] = [target]

        elif start in self.km1Map.keys():

            if target not in self.km1Map[start]:
                new = []
                for n in self.km1Map[start]:
                    new.append(n)
                new.append(target)
                self.km1Map[start] = new

class Bacteria():
    def __init__(self, name, id, match):
        self.name = name
        self.id = id
        self.match = match


    def __repr__(self):
        return f"<Bacteria {self.name}>"
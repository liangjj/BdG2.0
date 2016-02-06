from parameters import Parameters


class System:
    """The main interface with which to interact with the code.

        A system allows for calculating the superconducting gap, Tc, etc...
    """
    def __init__(self, material, shape):
        self.material = material
        self.shape = shape
        self.parameters = Parameters(material, shape)

    def get_tc():
        print("Hello")

#  This class contains various components used in electrical circuits. 
#  Component is a parent class for all the child "component" classes.
from math import acos
from Settings import settings


class Load:
    """
    Class to represent load objects
    """
    def __init__(self, name: str, bus: str, real_power: float, reactive_power: float):
        """
        Constructor for Load class
        :param name: Name of load
        :param bus: Bus connection
        :param real_power: Real power load is using
        :param reactive_power: Reactive power load is using
        """
        self.name = name
        self.bus = bus
        self.real_rated = real_power*1e6
        self.reactive_rated = reactive_power*1e6
        self.real_power = real_power*1e6
        self.reactive_power = reactive_power*1e6
        self.Smag = (self.real_power**2 + self.reactive_power**2)**(1/2)
        self.S = self.real_power + 1j*self.reactive_power
        self.pf = self.real_power/self.Smag
        self.angle = acos(self.pf)



class Generator:
    """
    Class to represent generator objects
    """
    def __init__(self, name: str, bus: str, voltage: float, real_power: float, sub_transient_reactance = 0.0, neg_impedance = 0.0, zero_impedance = 0.0, gnd_impedance = None, var_limit = float('inf')):
        """
        Constructor for Generator class
        :param name: Name of generator
        :param bus: Bus connection
        :param voltage: Generator voltage
        :param real_power: Real power generation
        :param sub_transient_reactance: Sub transient reactance for faults
        :param neg_impedance: Negative impedance
        :param zero_impedance: Zero impedance
        :param gnd_impedance: Ground impedance
        :param var_limit: VAR Limit
        """
        self.name = name
        self.bus = bus
        self.voltage = voltage
        self.real_power = real_power*1e6
        self.reactive_power = 0.
        self.X0 = self.calc_X0(zero_impedance)
        self.X1 = self.calc_X1(sub_transient_reactance)
        self.X2 = self.calc_X2(neg_impedance)
        self.Zn = gnd_impedance
        self.Y0prim = self.calc_Y0prim()
        self.var_limit = var_limit
    

    def calc_X0(self, X0):
        """
        Return imaginary reactance X0
        :param X0: Reactance
        :return:
        """
        X0 = 1j*X0*settings.powerbase/self.real_power
        return X0


    def calc_X1(self, X1):
        """
        Return imaginary reactance X1
        :param X1: Reactance
        :return:
        """
        X1 = 1j*X1*settings.powerbase/self.real_power
        return X1


    def calc_X2(self, X2):
        """
        Return imaginary reactance X2
        :param X2: Reactance
        :return:
        """
        X2 = 1j*X2*settings.powerbase/self.real_power
        return X2


    def set_power(self, real: float, reactive: float):
        """
        Set function for power
        :param real: Real power
        :param reactive: Reactive power
        :return:
        """
        self.real_power = real*1e6
        self.reactive_power = reactive*1e6

        
    def calc_Y0prim(self):
        """
        Generator primitive admittance matrix
        :return:
        """
        if self.Zn == None:
            Y0prim = 0
        elif self.Zn >= 0:
            Y0prim = 1/(3*self.Zn+self.X0)
        
        return Y0prim

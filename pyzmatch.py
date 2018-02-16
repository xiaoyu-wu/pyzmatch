# Standard Libraries
import math
import cmath
import numbers

# Major Libraries
import matplotlib.pyplot as plt
import numpy as np

########Global Viarables########
FREQ = 1  # frequency in GHz
################################

########Class Def###############
class Cable:

    def __init__(self, loss, velocity, length, impedance=50):
        # impedance in Ohm, loss in dB/100ft, velocity in c, length in cm
        # self.loss is transformed into /m
        self.impedance = impedance
        self.loss = - math.log (10 ** (- loss / 20.0)) / 30.5
        self.vel = velocity
        self.length = length / 100.0
        self.nexts = []

    def alpha(self):
        global FREQ
        alpha = self.loss * math.sqrt(FREQ)
        return alpha

    def beta(self):
        global FREQ
        beta = 2 * math.pi * FREQ * (10 ** 9) / (self.vel * 3 * (10 ** 8))
        return beta

    def gamma(self):
        gamma = complex(self.alpha(),self.beta())
        return gamma

    def set_length(self, length):
        self.length = length / 100

    def set_nexts(self, next_element):
        self.nexts.append(next_element)

    def network_impedance(self):
        if self.nexts == []:
            impedance = Z_cable(self)
        else:
            admittance = 0
            for i in self.nexts:
                admittance += 1.0/i.network_impedance()
            load = 1.0/admittance
            impedance = Z_distant(load, self)
        return impedance


class Lumped_Element:

    def __init__(self, resistance, capacitance, inductance):
        # PARALLEL resistance in Ohm, Capacitance in pF, inductance in nH
        self.resistance = resistance
        self.capacitance = capacitance
        self.inductance = inductance
        self.nexts = []

    def self_impedance(self):
        global FREQ
        if self.capacitance != 0:
            impedance = complex(self.resistance, 2 * math.pi * FREQ * self.inductance - 1.0 / ( 2 * math.pi * FREQ * self.capacitance * 10 ** (-3)))
        else:
            impedance = complex(self.resistance, 2 * math.pi * FREQ * self.inductance)
        return impedance

    def set_r(self, resistance):
        self.resistance = resistance

    def set_c(self, capacitance):
        self.capacitance = capacitance

    def set_l(self, inductance):
        self.inductance = inductance

    def set_nexts(self, next_element):
        self.nexts.append(next_element)

    def network_impedance(self):
        if (self.nexts == []):
            impedance = self.self_impedance()
        else:
            impedance = 0
            admittance = 0
            for i in self.nexts:
                admittance += 1.0/i.network_impedance()
            impedance = 1.0 / admittance + self.self_impedance()
        return impedance


###############################



#######Function Def############
def V_reflection(impedance1, impedance2):
    #voltage reflection when wave goes from media1 with impedance1 to media2 with impedance2
    tau = complex(impedance2 - impedance1) / complex(impedance2 + impedance1)
    return tau


def Z_distant(Z_load, cable):
    #calculate the impedance seen at one end of a cable with the other end connecting to a load with impedance Z_load
    Z_distant = cable.impedance * (1 + V_reflection(cable.impedance, Z_load) * cmath.exp(-2 * cable.gamma() * cable.length)) / (1 - V_reflection(cable.impedance, Z_load) * cmath.exp(-2 * cable.gamma() * cable.length))
    return Z_distant


def Z_cable(cable):
    #calculate the impedance seen at one end of a cable with the other end open
    #the same effect with Z_distant(10 ** 10, cable)
    Z_cable = cable.impedance * (1 + cmath.exp(-2 * cable.gamma() * cable.length)) / (1 - cmath.exp(-2 * cable.gamma() * cable.length))
    return Z_cable


def change_parameter(target, parameter_name, new_value):
    global FREQ
    if (parameter_name == "frequency"):
        FREQ = new_value
    else:
        dict_short = {"resistance":"r", "capacitance":"c", "inductance":"l", "length":"length"}
        getattr(target, "set_"+dict_short[parameter_name])(new_value)


def one_parameter_search(target, parameter_name, parameter_range, END, tip):
    S11_list = []
    dS11_list = []
    S11_linear = []
    f = open("./one_parameter_search_output.txt","w")

    for i in parameter_range:
        f.write(str(i)+",")
        change_parameter(target, parameter_name, i)
        #calculate S11
        S11 = V_reflection(50, END.network_impedance())
        S11_linear.append(abs(S11))
        f.write(str(abs(S11))+"\n")
        S11_dB = 20 * math.log(abs(S11), 10)
        S11_list.append(S11_dB)
        #calculate delta S11
        tip.set_c(tip.capacitance+10**-6)
        dS11 = abs(V_reflection(50, END.network_impedance()) - S11)
        dS11_list.append(dS11)
        tip.set_c(tip.capacitance-10**-6)
    f.close()

    #Plot the simulation results
    fig, ax1 = plt.subplots()
    ax1.plot(parameter_range, S11_list, 'b-')
    ax1.set_xlabel(parameter_name)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('S11', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')


    ax2 = ax1.twinx()
    ax2.plot(parameter_range, dS11_list, 'r-')
    ax2.set_ylabel('Delta S11', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    plt.show()
    return S11_list, dS11_list


def two_parameter_search(target1, parameter_name1, parameter_range1, target2, parameter_name2, parameter_range2, END, tip):
    S11_map = []
    dS11_map = []
    for i in parameter_range1:
        change_parameter(target1, parameter_name1, i)
        S11_line = []
        dS11_line = []
        for j in parameter_range2:
            change_parameter(target2, parameter_name2, j)
            #S11
            S11 = V_reflection(50, END.network_impedance())
            S11_dB = 20 * math.log(abs(S11), 10)
            S11_line.append(S11_dB)
            #delta S11
            tip.set_c(tip.capacitance+10**-6)
            dS11 = abs(V_reflection(50, END.network_impedance()) - S11)
            dS11_line.append(dS11)
            tip.set_c(tip.capacitance-10**-6)
        S11_map.append(S11_line)
        dS11_map.append(dS11_line)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    f, (ax1, ax2) = plt.subplots(ncols=2)
    im1 = ax1.imshow(S11_map, aspect='auto', extent=[parameter_range2[0], parameter_range2[-1], parameter_range1[0], parameter_range1[-1]],  origin='lower')
    ax1.set_title('S11')
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="3%", pad=0.1)
    cbar1 = plt.colorbar(im1, cax=cax1)
    im2 = ax2.imshow(dS11_map, aspect='auto', extent=[parameter_range2[0], parameter_range2[-1], parameter_range1[0], parameter_range1[-1]],  origin='lower')
    ax2.set_title('dS11')
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="3%", pad=0.1)
    cbar2 = plt.colorbar(im2, cax=cax2)
    plt.show()

###############################

##########Example##############

def test_stub_tuning():
    global FREQ
    #Define Z-Match Config.
    END = Lumped_Element(0, 0, 0)
    Stub = Cable(20, 0.7, 4.5)
    QWC = Cable(73, 0.7, 4.2)
    tip = Lumped_Element(4, 1, 2)
    END.set_nexts(Stub)
    END.set_nexts(QWC)
    QWC.set_nexts(tip)
    #Search!
    two_parameter_search(FREQ, "frequency", np.linspace(0.7, 1.1, 100), Stub, "length", np.linspace(4, 7, 100), END, tip)
    Stub.set_length(4.57)
    one_parameter_search(FREQ, "frequency", np.linspace(0.9, 1.1, 200), END, tip)

if __name__ == '__main__':
    test_stub_tuning()

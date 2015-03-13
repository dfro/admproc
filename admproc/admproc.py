# -*- coding: utf-8 -*-
import numpy


class Data(object):
    """Class for adm data file

     attributes:
     """
    def __init__(self, fname):
        self.data, self.freq, self.area, self.eps = read(fname)
        self.fsel = None
        self.vsel = None
        self.tsel = None
        if self.freq.shape[0] == 1:
            self.fsel = self.freq[0]

    def extract(self, fsel=None, vsel=None, tsel=None):
        if fsel is not None:
            self.fsel = fsel
        if vsel is not None:
            self.vsel = vsel
        if tsel is not None:
            self.tsel = tsel
        extr_data = extract(self.data, self.freq, self.fsel, self.vsel, self.tsel)
        self.cap = extr_data[0]
        self.cond = extr_data[1]
        self.voltage = extr_data[2]
        self.temp = extr_data[3]
        # dissipation factor
        self.diss = self.cond/(2*numpy.pi*self.fsel*self.cap)

    def cp(self, fsel=None, vsel=None, tsel=None):
        """exact parallel capacitance"""
        self.extract(fsel, vsel, tsel)
        return self.cap

    def cs(self, fsel=None, vsel=None, tsel=None):
        """exact serial capacitance"""
        self.extract(fsel, vsel, tsel)
        return self.cap*(1 + self.diss**2)


def _adm_header(file):
    """find column headers in admittance file"""
    area = 0
    eps = 0
    for i, line in enumerate(file, 1):
        if '#Epsilon' in line:
            # used in old version of admfile
            area = float(line.split('\t')[1].split('=')[1])
            eps = float(line.split('\t')[2].split('=')[1])
        if line.startswith('area ='):
            area = float(line.split('=')[1])
        if line.startswith('epsilon ='):
            eps = float(line.split('=')[1])
        if line.startswith('#Temp'):
            # read frequencies form header line
            freq = []
            for column in line.split('\t'):
                if column.startswith('#C'):
                    freq.append(int(''.join(ele for ele in column if ele.isdigit())))
            freq = numpy.array(freq)
            break
    else:
        raise ValueError('Header line not found. End of file reached')
    return i, freq, area, eps


def read(fname):
    """read data array from admittance file
    
    Parameters
    ----------
    fname : file or str, file name
    
    
    Returns
    -------
    data : numpy array
    freq : float, frequency
    area : float, area of contact
    eps : float, epsilon - relative permittivity of material
    
    """
    
    # Support of in memory archive extraction
    try:
        with open(fname) as file:
            header_line, freq, area, eps = _adm_header(file)
   
    except TypeError:
        _, freq, area, eps = _adm_header(fname)
        header_line = 0
    
    data = numpy.genfromtxt(fname, skip_header=header_line)
    return data, freq, area, eps


def extract(data, freq, fsel=None, vsel=None, tsel=None):
    """extract capacitance and conductance form data array. 
    NOTE: at least two of sel parameters should be specified.
    
    Parameters
    ----------
    fsel : float, selected frequency
    vsel : float, selected voltage
    tsel : float, selected temperature
    data : array from read_data function
    freq : array from read_data function
    
    Returns
    -------
    cap : float, capacitance
    cond : float, conductance
    voltage : float, voltage
    temp : float, temperature
    
    """
    cap = []
    cond = []

    # read voltages
    voltage = []
    single_temp = False
    for i in range(len(data[:, 0])):
        if data[i, 0] != data[0, 0]:
            break
        voltage.append(round(data[i, 1], 6))

    else:
        # only one temperature in file
        single_temp = True
        tsel = data[0, 0]
    voltage = numpy.array(voltage)
    
    # read temperature
    temp = data[0::voltage.shape[0], 0]

    if freq.shape[0] == 1:
        fsel = freq[0]
    if voltage.shape[0] == 1:
        vsel = voltage[0]

    if tsel is None:
        # read C-T and G-T
        i = numpy.argmin(abs(voltage - vsel))
        j = numpy.argmin(abs(freq - fsel))
        if not vsel and not fsel:
            raise ValueError('at least two parameters should be specified')
        cap = data[i::voltage.shape[0], j+2]
        cond = data[i::voltage.shape[0], j+freq.shape[0]+2]
    if vsel is None:
        # read C-V and G-V
        j = numpy.argmin(abs(freq - fsel))
        k = numpy.argmin(abs(temp - tsel))*voltage.shape[0]
        cap = data[k:k+voltage.shape[0], j+2]
        cond = data[k:k+voltage.shape[0], j+freq.shape[0]+2]
    if fsel is None:
        # read C-f and G-f
        i = numpy.argmin(abs(voltage - vsel))
        k = numpy.argmin(abs(temp - tsel))*voltage.shape[0]
        cap = data[k+i, 2:freq.shape[0]+2]
        cond = data[k+i, 2+freq.shape[0]:]
    return cap, cond, voltage, temp


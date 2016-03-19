from unittest import TestCase
import numpy as np
import admproc


class TestAdmproc(TestCase):
    def test_one_row(self):
        data, freq, area, eps = admproc.read('one_row.dat')
        freq_ref = np.array([20, 100, 1000])
        data_ref = [3E+2, -2E-1, 1.76841E-8, 7.52769E-9, 2.46918E-9,
                                 7.27778E-10, 2.08698E-10, 5.05154E-11]

        self.assertTrue(np.allclose(area, 0.1))
        self.assertTrue(np.allclose(eps, 10))
        self.assertTrue(np.allclose(freq, freq_ref))
        self.assertTrue(np.allclose(data, data_ref))

    def test_old(self):
        file = 'old_fileformat.dat'
        data, freq, area, eps = admproc.read(file)
        cap, cond, voltage, temp = admproc.extract(data, freq, fsel=100, tsel=6)

        freq_ref = np.array([1000000, 500000, 200000, 100000, 50000, 20000,
                             10000, 5000, 2000, 1000])
        voltage_ref = np.linspace(2.5, 0, 26)

        self.assertTrue(np.allclose(area, 0.00144))
        self.assertTrue(np.allclose(eps, 8.90))
        self.assertTrue(np.allclose(freq, freq_ref))
        self.assertTrue(np.allclose(voltage, voltage_ref))

    def test_one_temperature(self):
        file = 'one_temperature.dat'
        data, freq, area, eps = admproc.read(file)
        cap, cond, voltage, temp = admproc.extract(data, freq, fsel=100)

        freq_ref = np.array([20, 100, 1000])
        voltage_ref = np.linspace(-0.2, 2.2, 13)

        self.assertTrue(np.allclose(area, 0.1))
        self.assertTrue(np.allclose(eps, 10))
        self.assertTrue(np.allclose(freq, freq_ref))
        self.assertTrue(np.allclose(voltage, voltage_ref))
        self.assertTrue(np.allclose(temp, 300))

    def test_nxcalc(self):
        volt = np.array([0.75, 0.65, 0.55, 0.45, 0.35, 0.25])
        cap = np.array([0.11543373, 0.1259744, 0.14012403,
                        0.16048704, 0.19359813, 0.26199768])
        dop_ref = -1e17
        width_ref =np.array([85.82477936, 77.861412, 68.92227719,
                            58.5136042, 45.47627313])

        dop, width = admproc.nxcalc(cap*1e-6, volt, area=1, eps=11.7)

        self.assertTrue(np.allclose(width, width_ref))
        self.assertTrue(np.allclose(dop, dop_ref*np.ones(5),rtol=0.1))

    def test_openzip(self):
        pass

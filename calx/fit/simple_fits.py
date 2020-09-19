import numpy as np
from scipy.stats import skew
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

def peak_ana(x, y, nb=3, plotpoints=False):
    """ nb = number of point (on each side) to use as background"""
    ## get background
    xb = np.hstack((x[0:nb], x[-(nb):]))
    yb = np.hstack((y[0:nb], y[-(nb):]))
    a = np.polyfit(xb, yb, 1)
    b = np.polyval(a, x)
    yf = y - b
    yd = np.diff(yf)

    ## determine whether peak or step
    ispeak = np.abs(skew(yf)) > np.abs(skew(yd))
    if ispeak:
        yw = yf
        xw = x
    else:
        yw = yd
        xw = (x[1:] + x[0:-1]) / 2
        ## get background
        xwb = np.hstack((xw[0:nb], xw[-(nb):]))
        ywb = np.hstack((yw[0:nb], yw[-(nb):]))
        aw = np.polyfit(xwb, ywb, 1)
        bw = np.polyval(aw, xw)
        yw = yw - bw

    Iw = (xw[1:] - xw[0:-1]) * (yw[1:] + yw[0:-1]) / 2
    if sum(Iw) < 0:
        yw = -yw

    ## get parameters
    mm = yw.argmax(0)
    PEAK = xw[mm]
    ywmax = yw[mm]
    gg = (yw[:mm][::-1] < (ywmax / 2)).argmax()
    ip = interp1d(
        yw.take([mm - gg - 1, mm - gg]), xw.take([mm - gg - 1, mm - gg]), kind="linear"
    )
    xhm1 = ip(ywmax / 2)
    gg = (yw[mm:] < (ywmax / 2)).argmax()
    ip = interp1d(
        yw.take([mm + gg, mm + gg - 1]), xw.take([mm + gg, mm + gg - 1]), kind="linear"
    )
    xhm2 = ip(ywmax / 2)

    FWHM = np.abs(xhm2 - xhm1)
    CEN = (xhm2 + xhm1) / 2
    if plotpoints and ispeak is True:
        # plot the found points for center and FWHM edges
        ion()
        plt.hold(True)
        plt.plot(x, b, "g--")
        plt.plot(x, b + ywmax, "g--")
        plt.plot([xhm1, xhm1], polyval(a, xhm1) + [0, ywmax], "g--")
        plt.plot([xhm2, xhm2], polyval(a, xhm2) + [0, ywmax], "g--")
        plt.plot([CEN, CEN], polyval(a, CEN) + [0, ywmax], "g--")
        plt.plot([xhm1, xhm2], [polyval(a, xhm1), polyval(a, xhm2)] + ywmax / 2, "gx")
        plt.draw()

    if not ispeak:
        try:
            # findings start of step coming from left.
            std0 = sp.std(y[0:nb])
            nt = nb
            while (sp.std(y[0:nt]) < (2 * std0)) and (nt < len(y)):
                nt = nt + 1
            lev0 = sp.mean(y[0:nt])
            # findings start of step coming from right.
            std0 = sp.std(y[-nb:])
            nt = nb
            while (sp.std(y[-nt:]) < (2 * std0)) and (nt < len(y)):
                nt = nt + 1
            lev1 = sp.mean(y[-nt:])
            gg = np.abs(y - ((lev0 + lev1) / 2)).argmin()
            ftx = y[gg - 2 : gg + 2]
            fty = x[gg - 2 : gg + 2]
            if ftx[-1] < ftx[0]:
                ftx = ftx[::-1]
                fty = fty[::-1]
            ip = interp1d(ftx, fty, kind="linear")
            CEN = ip((lev0 + lev1) / 2)
            gg = np.abs(y - (lev1 + (lev0 - lev1) * 0.1195)).argmin()
            ftx = y[gg - 2 : gg + 2]
            fty = x[gg - 2 : gg + 2]
            if ftx[-1] < ftx[0]:
                ftx = ftx[::-1]
                fty = fty[::-1]
            # print " %f %f %f %f %f" % (ftx[0],ftx[1],fty[0],fty[1],lev1+(lev0-lev1)*0.1195)
            ip = interp1d(ftx, fty, kind="linear")
            H1 = ip((lev1 + (lev0 - lev1) * 0.1195))
            # print "H1=%f" % H1

            gg = np.abs(y - (lev0 + (lev1 - lev0) * 0.1195)).argmin()

            ftx = y[gg - 2 : gg + 2]
            fty = x[gg - 2 : gg + 2]

            if ftx[-1] < ftx[0]:
                ftx = ftx[::-1]
                fty = fty[::-1]
            #    print " %f %f %f %f %f" % (ftx[0],ftx[1],fty[0],fty[1],lev0+(lev1-lev0)*0.1195)
            ip = interp1d(ftx, fty, kind="linear")
            H2 = ip((lev0 + (lev1 - lev0) * 0.1195))
            # print "H2=%f" % abs(H2-H1)
            FWHM = abs(H2 - H1)
            if plotpoints is True:
                # plot the found points for center and FWHM edges
                plt.ion()
                plt.hold(True)
                plt.plot([x.min(), x.max()], [lev0, lev0], "g--")
                plt.plot([x.min(), x.max()], [lev1, lev1], "g--")
                plt.plot([H2, H2], [lev0, lev1], "g--")
                plt.plot([H1, H1], [lev0, lev1], "g--")
                plt.plot([CEN, CEN], [lev0, lev1], "g--")
                plt.plot(
                    [H2, CEN, H1],
                    [
                        lev0 + (lev1 - lev0) * 0.1195,
                        (lev1 + lev0) / 2,
                        lev1 + (lev0 - lev1) * 0.1195,
                    ],
                    "gx",
                )
                plt.draw()
        except:
            CEN = np.nan
            FWHM = np.nan
            PEAK = np.nan
    return (CEN, FWHM, PEAK)

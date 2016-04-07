package Modules;

import CrystTools.Fragment;
import CrystTools.Symmetry;
import CrystTools.UnitCell;
import MathTools.FastMath;
import Utilities.ObjectsUtilities;

import java.io.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.logging.ConsoleHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static MathTools.FastMath.*;
import static Utilities.ObjectsUtilities.*;

/**
 * Created by Developer on 16.09.2015.
 */
public class MoRmodule {

    public final UnitCell CELL;
    public FragmentData FRAG;
    public DiffractionData HKL;
    public Symmetry SYM;
    public RavinessSettings MORSETTINGS;
    private int[] ParametersList;
    private Energy E;
    private String tempFileName;
    private String resultFileName;
    private static Logger log = Logger.getLogger(MoRmodule.class.getName());

//    private static volatile boolean cancelled = false;
//    private Thread threadKeyboard = new Thread(new Runnable(){
//        public void run() {
//            Scanner scan = new Scanner(System.in);
//            while (true) {
//                if (scan.nextLine().equals("s")) cancelled = true;
//            }
//        }
//    });


    public class RavinessSettings {

        private String ravinessSettingsFilename;
        private int N = 0;
        private int REFRESH = 0;
        private double h = 2E-1;
        private double c1 = 1.0E-3;
        private double c2 = 1.0E-3;
        private double c3 = 1.0E-3;
        private double Epsilon1 = 1.0E-3;
        private double Epsilon2 = 1.0E-3;
        private double Delta = 0;
        private double Delta1 = 0;
        private double Delta2 = 0;
        private int N_TMIN = 100;
        private int FIRST_RANDOMIZATION = 50;
        private int MINIMA_SEARCH_METHOD = 1;
        private int N_FIB = 10;
        private String SAVE_BEST_RESULT;
        private int G_FIND_DERIVATIVE = 3;
        private int PRINT_FRAGMENT;
        private String Eopt = "";
        private String wExray;
        private String wEcore;
        private String K;

        //----------------------------------------------------------------------------------------------------------------------
        private RavinessSettings(String ravinessDataFilename) {
            this.ravinessSettingsFilename = ravinessDataFilename;
            List<String> input = ObjectsUtilities.getContentFromFile(this.ravinessSettingsFilename, "METHOD_OF_RAVINESS_SETTINGS");
            for (String s : input) {
                try {
                    Pattern p = Pattern.compile("(\\S+)");
                    Matcher m = p.matcher(s);
                    List<String> allMatches = new ArrayList<String>();
                    while (m.find()) allMatches.add(m.group());
                    switch (allMatches.get(0)) {
                        case "N":
                            this.N = Integer.parseInt(allMatches.get(1));
                            break;
                        case "REFRESH":
                            this.REFRESH = Integer.parseInt(allMatches.get(1));
                            break;
                        case "N_FIB":
                            int tmp = Integer.parseInt(allMatches.get(1));
                            this.N_FIB = (tmp > 99) ? 99 : tmp;
                            break;
                        case "h":
                            this.h = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "c1":
                            this.c1 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "c2":
                            this.c2 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "c3":
                            this.c3 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "Epsilon2":
                            this.Epsilon2 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "N_TMIN":
                            this.N_TMIN = Integer.parseInt(allMatches.get(1));
                            break;
                        case "Epsilon1":
                            this.Epsilon1 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "Delta2":
                            this.Delta2 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "Delta":
                            this.Delta = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "E":
                            this.Eopt = allMatches.get(1);
                            break;
                        case "wExray":
                            this.wExray = allMatches.get(1);
                            break;
                        case "wEcore":
                            this.wEcore = allMatches.get(1);
                            break;
                        case "K":
                            this.K = allMatches.get(1);
                            break;
                        case "Delta1":
                            this.Delta1 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "FIRST_RANDOMIZATION":
                            this.FIRST_RANDOMIZATION = Integer.parseInt(allMatches.get(1));
                            break;
                        case "G_FIND_DERIVATIVE":
                            this.G_FIND_DERIVATIVE = Integer.parseInt(allMatches.get(1));
                            break;
                        case "MINIMA_SEARCH_METHOD":
                            this.MINIMA_SEARCH_METHOD = Integer.parseInt(allMatches.get(1));
                            break;
                        case "SAVE_BEST_RESULT":
                            this.SAVE_BEST_RESULT = "";
                            break;
                        case "PRINT_FRAGMENT":
                            this.PRINT_FRAGMENT = Integer.parseInt(allMatches.get(1));
                            break;
                        default:
                            break;
                    }
                    allMatches.clear();
                } catch (NumberFormatException e) {
                    //throw new RuntimeException(e);
                }
            }
        }
    }
//----------------------------------------------------------------------------------------------------------------------


    public MoRmodule(UnitCell CELL,
                     Symmetry SYM,
                     FragmentData FRAG,
                     DiffractionData HKL,
                     Energy E,
                     String ravinessSettingsFilename,
                     String tempFileName,
                     String resultFileName) {
        this.CELL = CELL;
        this.SYM = SYM;
        this.FRAG = FRAG;
        this.HKL = HKL;
        this.MORSETTINGS = new RavinessSettings(ravinessSettingsFilename);
        this.tempFileName = tempFileName;
        this.resultFileName = resultFileName;
        DeleteFile(tempFileName);
        this.ParametersList = generateParametersList(FRAG);
        this.E = E;
        Handler consoleHandler = new ConsoleHandler();
        consoleHandler.setLevel(Level.FINER);
        log.getAnonymousLogger().addHandler(consoleHandler);
    }


    public double findt(FragmentData FRAG, Energy E, double[] g) {
        double t = this.MORSETTINGS.Delta1;
        int counter = this.MORSETTINGS.N_TMIN;
        E.calcEnergy(FRAG);
        double E0t = E.E;
        double E1t = 0;
        double E2t = 0;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-t, g));
        E.calcEnergy(FRAG);
        double mEt = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(2.0 * t, g));
        E.calcEnergy(FRAG);
        double pEt = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-t, g));
        do {
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(t, g));
            E.calcEnergy(FRAG);
            E1t = E.E;
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-t, g));
            if (E1t <= E0t) {
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(2.0 * t, g));
                E.calcEnergy(FRAG);
                E2t = E.E;
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-2.0 * t, g));
                if (E2t - 2 * E1t + E0t < 0) {
                    t *= 2.0;
                } else {
                    break;
                }
            } else {
                t /= 2.0;
            }
            counter--;
            if (counter == 0) {
                return this.MORSETTINGS.Delta1;
            }
        } while (true);
        return t;
    }


    public double findt2(FragmentData FRAG, Energy E, double[] g) {
        double t = this.MORSETTINGS.Delta1;
        while (findD2t(FRAG, E, g, t) >= 0) t *= 2;
        return t;
    }

//    public double findDt(FragmentData FRAG, Energy E, double[] g, double t) {
//        double dt = this.MORSETTINGS.DERIVATIVE_FIRST_STEP;
//        addToParameters(FRAG, FastMath.KmV(-dt, g));
//        E.calcEnergy(FRAG);
//        double E0 = E.E;
//        addToParameters(FRAG, FastMath.KmV(2 * dt, g));
//        E.calcEnergy(FRAG);
//        double E1 = E.E;
//        addToParameters(FRAG, FastMath.KmV(-dt, g));
//        double df = (E1 - E0) / dt / 2;
//        double df_pre = 0;
//        do {
//            df_pre = df;
//            dt /= 2;
//            addToParameters(FRAG, FastMath.KmV(-dt, g));
//            E.calcEnergy(FRAG);
//            E0 = E.E;
//            addToParameters(FRAG, FastMath.KmV(2 * dt, g));
//            E.calcEnergy(FRAG);
//            E1 = E.E;
//            addToParameters(FRAG, FastMath.KmV(-dt, g));
//            df = (E1 - E0) / dt / 2;
//        } while ((df_pre != df) && (dt > this.MORSETTINGS.c3));
//        return df;
//    }


    public double findDt(FragmentData FRAG, Energy E, double[] g, double t) {
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(t, g));
        double dt = this.MORSETTINGS.c3;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-dt, g));
        E.calcEnergy(FRAG);
        double E0 = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(2 * dt, g));
        E.calcEnergy(FRAG);
        double E1 = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-dt, g));
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-t, g));
        return (E1 - E0) / dt / 2.0;
    }


//    public double findD2t(FragmentData FRAG, Energy E, double[] g, double t) {
//        double dt = this.MORSETTINGS.DERIVATIVE_FIRST_STEP;
//        double d2f_pre = 0;
//        E.calcEnergy(FRAG);
//        double E1 = E.E;
//        addToParameters(FRAG, FastMath.KmV(-dt, g));
//        E.calcEnergy(FRAG);
//        double E0 = E.E;
//        addToParameters(FRAG, FastMath.KmV(2 * dt, g));
//        E.calcEnergy(FRAG);
//        double E2 = E.E;
//        addToParameters(FRAG, FastMath.KmV(-dt, g));
//        double d2f = (E0 - 2 * E1 + E2) / (dt * dt);
//        do {
//            d2f_pre = d2f;
//            dt /= 2;
//            addToParameters(FRAG, FastMath.KmV(-dt, g));
//            E.calcEnergy(FRAG);
//            E0 = E.E;
//            addToParameters(FRAG, FastMath.KmV(2 * dt, g));
//            E.calcEnergy(FRAG);
//            E2 = E.E;
//            addToParameters(FRAG, FastMath.KmV(-dt, g));
//            d2f = (E0 - 2 * E1 + E2) / (dt * dt);
//        } while ((d2f_pre != d2f) && (dt > this.MORSETTINGS.c3));
//        return d2f;
//    }


    public double findD2t(FragmentData FRAG, Energy E, double[] g, double t) {
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(t, g));
        double dt = this.MORSETTINGS.c3;
        E.calcEnergy(FRAG);
        double E1 = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-dt, g));
        E.calcEnergy(FRAG);
        double E0 = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(2 * dt, g));
        E.calcEnergy(FRAG);
        double E2 = E.E;
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-dt, g));
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-t, g));
        return (E0 - 2 * E1 + E2) / dt / dt;
    }

    public void circSort(double[] MASS) {
        double M = MASS[0];
        MASS[0] = MASS[1];
        MASS[1] = MASS[2];
        MASS[2] = M;
    }

    public void findMinLine(int mode, FragmentData FRAG, Energy E, double[] g) {
        double tmin = findt(FRAG, E, g);
        if (tmin == 0) return;
        if (mode == 1) {
            double Phi = (1 + Math.sqrt(5)) / 2.0;
            double l = 0;
            double r = tmin;
            double[] tV = {l + (r - l) / (Phi + 1), r - (r - l) / (Phi + 1)};
            double[] eV = {0, 0};
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[0], g));
            E.calcEnergy(FRAG);
            eV[0] = E.E;
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[1] - tV[0], g));
            E.calcEnergy(FRAG);
            eV[1] = E.E;
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-tV[1], g));
            do {
                if (eV[0] < eV[1]) {
                    r = tV[1];
                    tV[1] = tV[0];
                    eV[1] = eV[0];
                    tV[0] = l + (r - l) / (Phi + 1);
                    addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[0], g));
                    E.calcEnergy(FRAG);
                    eV[0] = E.E;
                    addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-tV[0], g));
                } else {
                    l = tV[0];
                    tV[0] = tV[1];
                    eV[0] = eV[1];
                    tV[1] = r - (r - l) / (Phi + 1);
                    addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[1], g));
                    E.calcEnergy(FRAG);
                    eV[1] = E.E;
                    addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-tV[1], g));
                }
            }
            while (Math.abs(r - l) > this.MORSETTINGS.Epsilon1);
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV((r + l) / 2, g));
        } else if (mode == 2) {
            double a = 0;
            double b = tmin;
            int N = this.MORSETTINGS.N_FIB;
            List<Double> F_R = FastMath.findFibonacciNumbersRatios(N);
            double x1 = 0;
            double x2 = 0;
            for (int i = 1; i < N - 1; i++) {
                x1 = a + F_R.get(N - i - 1 - 1) * (b - a);
                x2 = a + F_R.get(N - i - 1) * (b - a);
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(x1, g));
                E.calcEnergy(FRAG);
                double E1 = E.E;
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(x2 - x1, g));
                E.calcEnergy(FRAG);
                double E2 = E.E;
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-x2, g));
                if (E1 > E2) {
                    a = x1;
                } else if (E1 < E2) {
                    b = x2;
                } else if (E1 == E2) {
                    a = x1;
                    b = x2;
                }
                if (Math.abs(Math.abs(x1) - Math.abs(x2)) < this.MORSETTINGS.c3) break;
            }
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV((x1 + x2) / 2, g));

        } else if (mode == 3) {
            double t = tmin;
            double[] tV = new double[]{0, t, 2 * t};
            double[] eV = new double[]{0, 0, 0};
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[0], g));
            E.calcEnergy(FRAG);
            eV[0] = E.E;
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[1] - tV[0], g));
            E.calcEnergy(FRAG);
            eV[1] = E.E;
            addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[2] - tV[1], g));
            E.calcEnergy(FRAG);
            eV[2] = E.E;
            do {
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(-tV[2], g));
                double f01 = (eV[0] - eV[1]) / (tV[0] - tV[1]);
                double f02 = (eV[0] - eV[2]) / (tV[0] - tV[2]);
                double a = (f01 - f02) / (tV[1] - tV[2]);
                double b = f01 - a * (tV[0] + tV[1]);
                tV[0] = -b / (2 * a);
                circSort(tV);
                addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(tV[2], g));
                E.calcEnergy(FRAG);
                eV[0] = E.E;
                circSort(eV);
            }
            while (Math.abs(eV[1] - eV[2]) > this.MORSETTINGS.c1 * eV[2] | Math.abs(tV[1] - tV[2]) > this.MORSETTINGS.c2);
        }
    }


    public double[] findG(FragmentData FRAG, Energy E) {
        double[] dE = new double[this.ParametersList.length];
        double dx = this.MORSETTINGS.c3;
        for (int i = 0; i < this.ParametersList.length; i++) {
            dE[i] = partialDerivative(FRAG, E, dx, i, this.MORSETTINGS.G_FIND_DERIVATIVE);
        }
        double[] g = FastMath.KmV(-1.0, dE);
        double mdE = FastMath.modV(dE);
        for (int i = 0; i < this.ParametersList.length; i++) {
            g[i] = (mdE == 0) ? 0 : g[i] / mdE;
        }
        return g;
    }

    public double partialDerivative(FragmentData FRAG, Energy E, double dx, int i, int mode) {
        double result = 0.0;
        if (mode == 5) {
            addToParameters(CELL, FRAG, ParametersList, i, -2.0 * dx);
            E.calcEnergy(FRAG);
            double E0 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, +dx);
            E.calcEnergy(FRAG);
            double E1 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, +dx);
//            E.calcEnergy(FRAG);
//            double E2 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, +dx);
            E.calcEnergy(FRAG);
            double E3 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, +dx);
            E.calcEnergy(FRAG);
            double E4 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, -2.0 * dx);
            result = (E0 - 8.0 * E1 + 8.0 * E3 - E4) / dx / 12.0;
        }
        if (mode == 3) {
            addToParameters(CELL, FRAG, ParametersList, i, -dx);
            E.calcEnergy(FRAG);
            double E0 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, +2.0 * dx);
            E.calcEnergy(FRAG);
            double E2 = E.E;
            addToParameters(CELL, FRAG, ParametersList, i, -dx);
            result = (E2 - E0) / dx / 2.0;
        }
        return result;
    }


    public double[] findModsGEparts(FragmentData FRAG, Energy E) {
        double[] dEx = new double[this.ParametersList.length];
        double[] dEr = new double[this.ParametersList.length];
        double[] dEc = new double[this.ParametersList.length];
        double[] dEp = new double[this.ParametersList.length];
        double[] ModsG = new double[4];
        for (int i = 0; i < this.ParametersList.length; i++) {
            double dx = this.MORSETTINGS.c3;
            addToParameters(CELL, FRAG, ParametersList, i, -dx);
            E.calcEnergy(FRAG);
            double E0x = E.Exray;
            double E0r = E.Erest;
            double E0c = E.Ecore;
            double E0p = E.Epenalty;
            addToParameters(CELL, FRAG, ParametersList, i, 2 * dx);
            E.calcEnergy(FRAG);
            double E1x = E.Exray;
            double E1r = E.Erest;
            double E1c = E.Ecore;
            double E1p = E.Epenalty;
            addToParameters(CELL, FRAG, ParametersList, i, -dx);
            dEx[i] = (E1x - E0x) / dx / 2.0;
            dEr[i] = (E1r - E0r) / dx / 2.0;
            dEc[i] = (E1c - E0c) / dx / 2.0;
            dEp[i] = (E1p - E0p) / dx / 2.0;
        }
        ModsG[0] = Math.pow(FastMath.modV(dEx), 2);
        ModsG[1] = Math.pow(FastMath.modV(dEr), 2);
        ModsG[2] = Math.pow(FastMath.modV(dEc), 2);
        ModsG[3] = Math.pow(FastMath.modV(dEp), 2);
        return ModsG;
    }

    public double calcW(List<Double> E0, List<Double> E1, double w) {
        double E0Mid = 0;
        double E1Mid = 0;
        for (int i = 0; i < E0.size(); i++) E0Mid += E0.get(i);
        E0Mid /= E0.size();
        for (int i = 0; i < E1.size(); i++) E1Mid += E1.get(i);
        E1Mid /= E1.size();
        return (E0Mid * E1Mid != 0) ? 0.5 * Math.sqrt(E1Mid / E0Mid) : w;
    }


    public void gradientDec(FragmentData FRAG, Energy E, double Delta) {
        double E0 = 0;
        E.calcEnergy(FRAG);
        double E1 = E.E;
        do {
            E0 = E1;
            double[] g = findG(FRAG, E);
            if (FastMath.modV(g) == 0) break;
            findMinLine(this.MORSETTINGS.MINIMA_SEARCH_METHOD, FRAG, E, g);
            E.calcEnergy(FRAG);
            E1 = E.E;
        } while (E1 / E0 <= Delta);
    }

    public double[] findRho(FragmentData FRAGbegin, FragmentData FRAGend) {
        double[] rho = new double[this.ParametersList.length];
        for (int i = 0; i < this.ParametersList.length; i++) {
            rho[i] = getParameters(FRAGend, ParametersList, i) - getParameters(FRAGbegin, ParametersList, i);
        }

        return FastMath.KmV(1.0 / FastMath.modV(rho), rho);
    }

    public void jumpRaviness(FragmentData FRAG, double[] rho, double h) {
        addToParameters(CELL, FRAG, ParametersList, FastMath.KmV(h, rho));
    }

    public void printTempInfo(FragmentData FRAG, Energy E, int loop, String regime){
        List<String> output = new ArrayList<>();
        if(regime.equals("")) {
            output.add(String.format("%5s %4s %12s %6s %6s %6s %6s %7s %7s %7s %7s %7s %7s %7s %7s",
                    "#",
                    "TYPE",
                    "Ecore",
                    "RI",
                    "RII",
                    "RIII",
                    "RIV",
                    "X",
                    "Y",
                    "Z",
                    "Phi1",
                    "Phi2",
                    "Theta",
                    "O",
                    "U"));
        } else {
            Fragment frag = null;
            for (Fragment itemFrag : FRAG.getFragMass()) {
                if (itemFrag.getFragNum() == this.MORSETTINGS.PRINT_FRAGMENT) {
                    frag = itemFrag;
                    break;
                }
            }
            output.add(String.format("%5d %4s %8e %6.3f %6.3f %6.3f %6.3f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f",
                    loop,
                    regime,
                    E.Ecore,
                    E.RI,
                    E.RII,
                    E.RIII,
                    E.RIV,
                    frag.getFragX(),
                    frag.getFragY(),
                    frag.getFragZ(),
                    frag.getFragPhi1(),
                    frag.getFragPhi2(),
                    frag.getFragTheta(),
                    frag.getFragO(),
                    frag.getFragU()));
        }
        ObjectsUtilities.putContentToFile(this.tempFileName, output, true);
    }


    public void run() {

        List<Double> statE = new ArrayList<>();
        List<Double> statEcore = new ArrayList<>();
        List<Double> statModGsqExray = new ArrayList<>();
        List<Double> statModGsqErest = new ArrayList<>();
        List<Double> statK = new ArrayList<>();
        List<Double> statModGsqEpenalty = new ArrayList<>();
        List<Double> statModGsqEcore = new ArrayList<>();
        List<String> opt = new ArrayList<>();

        StringBuilder strOUT = new StringBuilder();
        StringBuilder strSTATUS = new StringBuilder();
        StringBuilder strIND = new StringBuilder();

        int numDecreasesGlobal = 0;
        double minEcore = 0;
        double minE = 0;


        strOUT.append(String.format("\n%-50s\n", "Method of raviness"));
        strOUT.append(String.format("%-50s\n", "+-----------------------------------------------------------------+"));
        strOUT.append(String.format("|  %-50s%12d |\n", "Number of fragments: ", FRAG.getFragMass().size()));
        strOUT.append(String.format("|  %-50s%12d |\n", "Number of used reflections: ", HKL.getHKL().size()));
        strOUT.append(String.format("|  %-50s%12d |\n", "Number of parameters: ", this.ParametersList.length));
        strOUT.append(String.format("|  %-50s%12.2e |\n", "Ravine jumping step: ", this.MORSETTINGS.h));
        strOUT.append(String.format("|  %-50s%12.2e |\n", "Ravine probe: ", this.MORSETTINGS.Delta));
        strOUT.append(String.format("|  %-50s%12s |\n", "Energy options: ", this.MORSETTINGS.Eopt));

        switch (this.MORSETTINGS.MINIMA_SEARCH_METHOD) {
            case 1:
                strOUT.append(String.format("|  %-50s%12s |\n", "Minima search method: ", "golden ratio"));
                break;
            case 2:
                strOUT.append(String.format("|  %-50s%12s |\n", "Minima search method: ",
                        "Fibonacci/" + String.valueOf(this.MORSETTINGS.N_FIB)));
                break;
            case 3:
                strOUT.append(String.format("|  %-50s%12s |\n", "Minima search method: ", "parabolic"));
                break;
            default:
                break;
        }
        strOUT.append(String.format("%-50s\n", "+-----------------------------------------------------------------+"));
        System.out.print(strOUT.toString());
        strOUT.setLength(0);

        if (this.ParametersList.length == 0) {
            E.OPT = MORSETTINGS.Eopt;
            E.calcEnergy(FRAG);
            E.printInfo();
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (S) - RI", E.RI));
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (C) - RII", E.RII));
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (W) - RIII", E.RIII));
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor - RIV", E.RIV));
            strOUT.append(String.format(" %-30s= %-12e\n", "Chem part of minimum E", E.Erest));
            strOUT.append(String.format(" %-30s= %-12e\n", "Energy", E.E));
            strOUT.append(String.format(" %-30s= %-12e\n", "Penalty function", E.Epenalty));
            System.out.print("\n" + strOUT + "\n");
            strOUT.setLength(0);
            strOUT.append(String.format("\n%-50s\n", "+-----------------------------------------------------------------+"));
            System.out.print(strOUT);
            System.out.print("Done.\n");
            FRAG.printFragsWithSym(SYM);
            System.exit(0);
        }


        System.out.println("Running...\n");
        if ((new File(this.resultFileName)).exists() && (this.MORSETTINGS.SAVE_BEST_RESULT != null)) {
            try {
                System.out.print("Fragment data loading...");
                FragmentData FRAG_LOAD = (FragmentData) recallObject(this.resultFileName);
                for (Fragment itemFrag : FRAG_LOAD.getFragMass()) {
                    itemFrag.setFragOPT(FRAG.getFragMass().get(FRAG_LOAD.getFragMass().indexOf(itemFrag)).getFragOPT());
                }
                FRAG = (FragmentData) deepClone(FRAG_LOAD);
                FRAG_LOAD = null;
            } catch (IOException | ClassNotFoundException e) {

            }
            System.out.print("\r\r");
        } else {
            randomizeParametersInitial(this.CELL, FRAG, MORSETTINGS.FIRST_RANDOMIZATION);
        }
        FRAG.printFrags();

        double[] rho;
        E.OPT = MORSETTINGS.Eopt;
        double h = MORSETTINGS.h;
        double Delta = MORSETTINGS.Delta;
        FragmentData FRAG_I = (FragmentData) deepClone(FRAG);
        FragmentData FRAG_II = (FragmentData) deepClone(FRAG);
        //FragmentData FRAG_III = (FragmentData) deepClone(FRAG);

        if (MORSETTINGS.PRINT_FRAGMENT != 0) printTempInfo(FRAG_I, E, 0, "");



        for (int i = 0; i < MORSETTINGS.N; i++) {

            //FRAG_I.fragsParametersAdjustment();

            if (this.MORSETTINGS.REFRESH - i % this.MORSETTINGS.REFRESH == 1) {
                statE.clear();
                statEcore.clear();
                statK.clear();
                statModGsqExray.clear();
                statModGsqErest.clear();
                statModGsqEcore.clear();
                statModGsqEpenalty.clear();
            }

            gradientDec(FRAG_I, E, Delta);
            E.calcEnergy(FRAG_I);
            if (this.MORSETTINGS.PRINT_FRAGMENT != 0) printTempInfo(FRAG_I, E, i + 1, "D");
            //E.calcEnergy(FRAG_I);
            statE.add(E.E);
            statEcore.add(E.Ecore);
            statK.add(E.getK());

            if (E.E < minE) {
                numDecreasesGlobal++;
                strIND.setLength(0);
                strIND.append(String.format(" %05d", numDecreasesGlobal));
                strIND.append(String.format(" RI=%-4.2f", E.RI));
                strIND.append(String.format(" RII=%-4.2f", E.RII));
                strIND.append(String.format(" RIII=%-4.2f", E.RIII));
                strIND.append(String.format(" RIV=%-4.2f", E.RIV));
                FRAG = (FragmentData) deepClone(FRAG_I);
                FRAG.printFrags();
                E.printInfo();
                if (!this.resultFileName.equals("") && (this.MORSETTINGS.SAVE_BEST_RESULT != null)) try {
                    saveObject(FRAG, this.resultFileName);
                } catch (IOException e) {
                }
            }
            //minEcore = statEcore.get(statEcore.indexOf(Collections.min(statEcore)));
            //minEcore = (minEcore == 0) ? statEcore.get(0) : Math.min(statEcore.get(i), minEcore);
            minEcore = (minEcore == 0) ? E.Ecore : Math.min(E.Ecore, minEcore);
            minE = (minE == 0) ? E.E : Math.min(E.E, minE);

            if ((MORSETTINGS.wExray != null) || (MORSETTINGS.wEcore != null)) {
                double[] ModsGparts = findModsGEparts(FRAG_I, E);
                statModGsqExray.add(ModsGparts[0]);
                statModGsqErest.add(ModsGparts[1]);
                statModGsqEcore.add(ModsGparts[2]);
                statModGsqEpenalty.add(ModsGparts[3]);
            }
            if (MORSETTINGS.wExray != null) {
                E.improvewExray(Double.valueOf(MORSETTINGS.wExray));
            } else {
                E.improvewExray(calcW(statModGsqExray, statModGsqErest, E.getwExray()));
            }
            if (MORSETTINGS.wEcore != null) {
                E.improvewEcore(Double.valueOf(MORSETTINGS.wEcore));
            } else {
                E.improvewEcore(calcW(statModGsqEcore, statModGsqEpenalty, E.getwEcore()));
            }

            if (MORSETTINGS.K != null) {
                E.improveK(Double.valueOf(MORSETTINGS.K));
            } else {
                E.improveK(FRAG_I);
            }

            rho = findRho(FRAG_II, FRAG_I);
            FRAG_II = (FragmentData) deepClone(FRAG_I);
            jumpRaviness(FRAG_I, rho, h);


            E.calcEnergy(FRAG_I);
            if (this.MORSETTINGS.PRINT_FRAGMENT != 0) printTempInfo(FRAG_I, E, i + 1, "J");

            strSTATUS.setLength(0);
            strSTATUS.append("Iteration number  ");
            System.out.print("\r" +
                    strSTATUS.substring(0, strSTATUS.length() - 2) +
                    ":" +
                    String.format(" %-5d ", i + 1) +
                    strIND);

        }
        strOUT.append(String.format("\n%-50s\n", "+-----------------------------------------------------------------+"));
        System.out.print(strOUT);
        System.out.print("Done.\n");

    }

}

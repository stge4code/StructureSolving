package Modules;


import CrystTools.*;
import MathTools.ComplexNumber;
import MathTools.FastMath;
import Utilities.ObjectsUtilities;
import org.w3c.dom.DOMImplementation;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by Developer on 11.09.2015.
 */
public class Energy {
    private DiffractionData HKL;
    private final UnitCell CELL;
    private final Symmetry SYM;
    public String OPT;
    public double E;
    public double Exray;
    public double Erest;
    public double Ecore;
    public double Epenalty;
    private double wExray;
    private double wEcore;
    private double K;
    public double RI;
    public double RII;
    public double RIII;
    public double RIV;
    private boolean autoAdjustK = false;
    private PenaltyFunction PSI;
    private String infoFileName;
    private EnergySettings ENERGYSETTINGS;

    private class EnergySettings {

        private String energySettingsFilename;
        private int WEIGHT_SCHEME = 0;
        private String PAIR_GENERATOR;
        private double A2 = 1.0;
        private double B2 = 1.0;
        private double K2 = 1.0;
        private double A4 = 1.0;
        private double B4 = 1.0;

        public EnergySettings(String energyDataFilename) {
            this.energySettingsFilename = energyDataFilename;
            List<String> input = ObjectsUtilities.getContentFromFile(this.energySettingsFilename, "ENERGY_SETTINGS");
            for (String s : input) {
                try {
                    Pattern p = Pattern.compile("(\\S+)");
                    Matcher m = p.matcher(s);
                    List<String> allMatches = new ArrayList<String>();
                    while (m.find()) allMatches.add(m.group());
                    switch (allMatches.get(0)) {
                        case "WEIGHT_SCHEME":
                            this.WEIGHT_SCHEME = Integer.parseInt(allMatches.get(1));
                            break;
                        case "A2":
                            this.A2 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "B2":
                            this.B2 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "A4":
                            this.A4 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "B4":
                            this.B4 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "K2":
                            this.K2 = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "PAIR_GENERATOR":
                            this.PAIR_GENERATOR = allMatches.get(1);
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


    public class Pair {
        public double Xi = 0;
        public double Yi = 0;
        public double Zi = 0;
        public double Xj = 0;
        public double Yj = 0;
        public double Zj = 0;
        public double R = 0;
        public double Rconstraint = 0;

        public Pair(double xi, double yi, double zi, double xj, double yj, double zj, double R, double Rconstraint) {
            Xi = xi;
            Yi = yi;
            Zi = zi;
            Xj = xj;
            Yj = yj;
            Zj = zj;
            this.R = R;
            this.Rconstraint = Rconstraint;
        }
    }


    public Energy(DiffractionData HKL,
                  UnitCell CELL,
                  Symmetry SYM,
                  PenaltyFunction PSI,
                  String energySettingsFilename,
                  String infoFilename) {
        this.ENERGYSETTINGS = new EnergySettings(energySettingsFilename);
        this.infoFileName = infoFilename;
        this.HKL = HKL;
        this.CELL = CELL;
        this.SYM = SYM;
        this.PSI = PSI;
        this.E = 0.0;
        this.Exray = 0.0;
        this.Erest = 0.0;
        this.Ecore = 0.0;
        this.Epenalty = 0.0;
        this.wExray = 1.0;
        this.wEcore = 1.0;
        this.K = 1.0;
        this.RI = 0.0;
        this.RII = 0.0;
        this.RIII = 0.0;
        this.RIV = 0.0;
    }


    public double getK() {
        return K;
    }

    public double getwExray() {
        return wExray;
    }

    public double getwEcore() {
        return wEcore;
    }

    public boolean ifInCell(double x, double y, double z) {
        boolean condition = (x <= 1) && (x >= 0) && (y <= 1) && (y >= 0) && (z <= 1) && (z >= 0);
        if (condition) return true;
        return false;
    }

    public double restrainFunction1(FragmentData FRAG) {
        List<Pair> Pairs = genAtomsPairs(this.SYM, FRAG);
        double restrain = 0;
//        double c = 1;
//        double n = 2;
//        double m = 2;
//              if ((R != 0) && (R < itemApair.RVdW))
//                  partErest += Math.pow(Math.pow(c * RVdW, n) - Math.pow(R, n), m);
        for (Pair itempair : Pairs) {
            if ((itempair.R != 0) && (itempair.R < itempair.Rconstraint)) {
                restrain += Math.pow(Math.min(itempair.R - itempair.Rconstraint, 0), 2);
            }
        }
        return restrain;
    }


    public double restrainFunction2(FragmentData FRAG) {
        List<Pair> Pairs = genFragmentsPairs(this.SYM, FRAG);
        double restrain = 0;
//        double c = 1;
//        double n = 2;
//        double m = 2;
//              if ((R != 0) && (R < itemApair.RVdW))
//                  partErest += Math.pow(Math.pow(c * RVdW, n) - Math.pow(R, n), m);
        for (Pair itempair : Pairs) {
            if ((itempair.R != 0) && (itempair.R < itempair.Rconstraint)) {
                restrain += Math.pow(Math.min(itempair.R - itempair.Rconstraint, 0), 2);
            }
        }
        return restrain;
    }


    public List<Pair> genAtomsPairs(Symmetry SYM, FragmentData FRAG) {
        List<Pair> result = new ArrayList<>();
        for (Fragment imFragment_i : FRAG.getFragMass()) {
            for (Fragment imFragment_j : FRAG.getFragMass()) {
                for (Atom imAtom_i : imFragment_i.getFragAtoms()) {
                    for (Atom imAtom_j : imFragment_j.getFragAtoms()) {
                        double Rconstraint = (imAtom_i.getAtomRVdW() + imAtom_j.getAtomRVdW());
                        for (SymmetryItem imSYM_i : SYM.getSymMass()) {
                            double[] VECT_i =
                                    FastMath.VpV(imSYM_i.getSymT(),
                                            FastMath.MmV(imSYM_i.getSymR(),
                                                    new double[]{
                                                            imAtom_i.getAtomX(),
                                                            imAtom_i.getAtomY(),
                                                            imAtom_i.getAtomZ()
                                                    }));
                            for (SymmetryItem imSYM_j : SYM.getSymMass()) {
                                double[] VECT_j =
                                        FastMath.VpV(imSYM_j.getSymT(),
                                                FastMath.MmV(imSYM_j.getSymR(),
                                                        new double[]{
                                                                imAtom_j.getAtomX(),
                                                                imAtom_j.getAtomY(),
                                                                imAtom_j.getAtomZ()
                                                        }));
                                if ((this.ENERGYSETTINGS.PAIR_GENERATOR != null) && this.ENERGYSETTINGS.PAIR_GENERATOR.contains("T")) {
                                    double[] VECT_it = VECT_i;
                                    for (int tX = -1; tX < 2; tX++) {
                                        for (int tY = -1; tY < 2; tY++) {
                                            for (int tZ = -1; tZ < 2; tZ++) {
                                                double[] VECT_jt = FastMath.VpV(VECT_j, new double[]{tX, tY, tZ});
                                                /*if (ifInCell(
                                                        (VECT_jt[0] + VECT_it[0]) / 2.0,
                                                        (VECT_jt[1] + VECT_it[1]) / 2.0,
                                                        (VECT_jt[2] + VECT_it[2]) / 2.0)) */{
                                                    double R = this.CELL.calcDistance(
                                                            VECT_it[0] - VECT_jt[0],
                                                            VECT_it[1] - VECT_jt[1],
                                                            VECT_it[2] - VECT_jt[2]);
                                                    result.add(new Pair(
                                                            VECT_it[0], VECT_it[1], VECT_it[2],
                                                            VECT_jt[0], VECT_jt[1], VECT_jt[2], R, Rconstraint));
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    /*if (ifInCell(
                                            (VECT_j[0] + VECT_i[0]) / 2.0,
                                            (VECT_j[1] + VECT_i[1]) / 2.0,
                                            (VECT_j[2] + VECT_i[2]) / 2.0))*/ {
                                        double R = this.CELL.calcDistance(
                                                VECT_i[0] - VECT_j[0],
                                                VECT_i[1] - VECT_j[1],
                                                VECT_i[2] - VECT_j[2]);
                                        result.add(new Pair(
                                                VECT_i[0], VECT_i[1], VECT_i[2],
                                                VECT_j[0], VECT_j[1], VECT_j[2], R, Rconstraint));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    }


    public List<Pair> genFragmentsPairs(Symmetry SYM, FragmentData FRAG) {
        List<Pair> result = new ArrayList<>();
        for (Fragment imFragment_i : FRAG.getFragMass()) {
            for (Fragment imFragment_j : FRAG.getFragMass()) {
                double Rconstraint = (imFragment_i.getFragDiameter() + imFragment_j.getFragDiameter()) / 2.0;
                for (SymmetryItem imSYM_i : SYM.getSymMass()) {
                    double[] VECT_i =
                            FastMath.VpV(imSYM_i.getSymT(),
                                    FastMath.MmV(imSYM_i.getSymR(),
                                            new double[]{
                                                    imFragment_i.getFragX(),
                                                    imFragment_i.getFragY(),
                                                    imFragment_i.getFragZ()
                                            }));
                    for (SymmetryItem imSYM_j : SYM.getSymMass()) {
                        double[] VECT_j =
                                FastMath.VpV(imSYM_j.getSymT(),
                                        FastMath.MmV(imSYM_j.getSymR(),
                                                new double[]{
                                                        imFragment_j.getFragX(),
                                                        imFragment_j.getFragY(),
                                                        imFragment_j.getFragZ()
                                                }));
                        if ((this.ENERGYSETTINGS.PAIR_GENERATOR != null) && this.ENERGYSETTINGS.PAIR_GENERATOR.contains("T")) {
                            double[] VECT_it = VECT_i;
                            for (int tX = -1; tX < 2; tX++) {
                                for (int tY = -1; tY < 2; tY++) {
                                    for (int tZ = -1; tZ < 2; tZ++) {
                                        double[] VECT_jt = FastMath.VpV(VECT_j, new double[]{tX, tY, tZ});
                                        /*if (ifInCell(
                                                (VECT_jt[0] + VECT_it[0]) / 2.0,
                                                (VECT_jt[1] + VECT_it[1]) / 2.0,
                                                (VECT_jt[2] + VECT_it[2]) / 2.0))*/ {
                                            double R = this.CELL.calcDistance(
                                                    VECT_it[0] - VECT_jt[0],
                                                    VECT_it[1] - VECT_jt[1],
                                                    VECT_it[2] - VECT_jt[2]);
                                            result.add(new Pair(
                                                    VECT_it[0], VECT_it[1], VECT_it[2],
                                                    VECT_jt[0], VECT_jt[1], VECT_jt[2], R, Rconstraint));
                                        }
                                    }
                                }
                            }
                        } else {
                            /*if (ifInCell(
                                    (VECT_j[0] + VECT_i[0]) / 2.0,
                                    (VECT_j[1] + VECT_i[1]) / 2.0,
                                    (VECT_j[2] + VECT_i[2]) / 2.0))*/ {
                                double R = this.CELL.calcDistance(
                                        VECT_i[0] - VECT_j[0],
                                        VECT_i[1] - VECT_j[1],
                                        VECT_i[2] - VECT_j[2]);
                                result.add(new Pair(
                                        VECT_i[0], VECT_i[1], VECT_i[2],
                                        VECT_j[0], VECT_j[1], VECT_j[2], R, Rconstraint));
                            }
                        }
                    }
                }
            }
        }
        return result;
    }


    public double calcwHKL(ReciprocalItem itemHKL) {
        double wHKL = 1.0;
        switch (this.ENERGYSETTINGS.WEIGHT_SCHEME) {
            case 0:
                wHKL = 1.0;
                break;
            case 1:
                double A1 = 2.0 * this.HKL.getParameters().IMin;
                double B1 = 2.0 / this.HKL.getParameters().IMax;
                wHKL = 1.0 / (A1 + itemHKL.Fo.getModule() + B1 * Math.pow(itemHKL.Fo.getModule(), 2));
                break;
            case 2:
                double A2 = this.ENERGYSETTINGS.A2;
                double B2 = this.ENERGYSETTINGS.B2;
                double K2 = this.ENERGYSETTINGS.K2;
                wHKL = 1.0 / (1.0 + Math.pow((K2 * itemHKL.Fo.getModule() - A2) / B2, 2));
                break;
            case 3:
                wHKL = 1.0 / Math.pow(itemHKL.Fo.getModule(), 2);
                break;
            case 4:
                double A4 = this.ENERGYSETTINGS.A4;
                double B4 = this.ENERGYSETTINGS.B4;
                if (!Double.isNaN(itemHKL.Fc.getModule())) {
                    double P = (2.0 * Math.pow(itemHKL.Fc.getModule(), 2) +
                            Math.max(Math.pow(itemHKL.Fo.getModule(), 2), 0.0)) / 3.0;
                    wHKL = 1.0 / (Math.pow(itemHKL.sigmaI, 2) + Math.pow(A4 * P, 2) + B4 * P);
                }
                break;
            default:
                break;
        }
        return wHKL;
    }


    public double PattersonFunction(double[] UVW) {
        double V0 = this.CELL.getV();
        double sumPatt = 0;
        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            sumPatt += itemHKL.I * Math.cos(2 * Math.PI * (itemHKL.h * UVW[0] + itemHKL.k * UVW[1] + itemHKL.l * UVW[2]));
        }
        return sumPatt / V0;
    }

    public double PattersonEnergy(Symmetry SYM, FragmentData FRAG) {
        double patt_E = 0;
        for (SymmetryItem imSYM : SYM.getSymMass()) {
            for (Fragment imFragment : FRAG.getFragMass()) {
                for (Atom imAtom_i : imFragment.getFragAtoms()) {
                    for (Atom imAtom_j : imFragment.getFragAtoms()) {
                        double ZiZj = imAtom_i.getAtomM() * imAtom_j.getAtomM();
                        double[] VECT_i =
                                FastMath.VpV(imSYM.getSymT(),
                                        FastMath.MmV(imSYM.getSymR(),
                                                new double[]{
                                                        imAtom_i.getAtomX(),
                                                        imAtom_i.getAtomY(),
                                                        imAtom_i.getAtomZ()
                                                }));
                        double[] VECT_j =
                                FastMath.VpV(imSYM.getSymT(),
                                        FastMath.MmV(imSYM.getSymR(),
                                                new double[]{
                                                        imAtom_j.getAtomX(),
                                                        imAtom_j.getAtomY(),
                                                        imAtom_j.getAtomZ()
                                                }));
                        double[] UVW =
                                {(VECT_j[0] - VECT_i[0]), (VECT_j[2] - VECT_i[1]), (VECT_j[2] - VECT_i[2])};
                        patt_E += Math.pow(ZiZj - PattersonFunction(UVW), 2);
                    }
                }
            }
        }
        return patt_E;
    }

    public void printInfo() {
        List<String> output = new ArrayList<>();
        output.add(String.format("%4s %4s %4s %14s %14s %14s %14s %14s %14s",
                "h",
                "k",
                "l",
                "SCATVECT",
                "Mod(Fo)",
                "Mod(Fc)*K",
                "Mod(Fc)",
                "Arg(Fc)",
                "Io-K*Ic"));
        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            output.add(String.format("%4d %4d %4d % .7e % .7e % .7e % .7e % .7e % .7e",
                    itemHKL.h,
                    itemHKL.k,
                    itemHKL.l,
                    FastMath.round(itemHKL.scatvect, 8),
                    FastMath.round(itemHKL.Fo.getModule(), 8),
                    FastMath.round(itemHKL.Fc.getModule() * this.K, 8),
                    FastMath.round(itemHKL.Fc.getModule(), 8),
                    FastMath.round(itemHKL.Fc.getPhase(), 8),
                    FastMath.round(Math.pow(itemHKL.Fo.getModule(), 2) - Math.pow(this.K * itemHKL.Fc.getModule(), 2), 8)));
        }
        ObjectsUtilities.putContentToFile(this.infoFileName, output);
    }

    public void setAutoImprovementK(boolean autoAdjustK) {
        this.autoAdjustK = autoAdjustK;
    }

    public void improveK(FragmentData FRAG) {
        double sumFoFc = 0;
        double sumFcFc = 0;
        double sumFcIm = 0;
        double sumFcRe = 0;

        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            double wHKL = calcwHKL(itemHKL);
            sumFcIm = sumFcRe = 0;
            for (Fragment itemFrag : FRAG.getFragMass()) {
                ComplexNumber fScat = itemFrag.fragScattering(itemHKL, this.CELL, this.SYM);
                sumFcRe += fScat.getRe();
                sumFcIm += fScat.getIm();
            }
            itemHKL.Fc.setNum(sumFcRe, sumFcIm);
            sumFoFc += wHKL * itemHKL.Fo.getModule() * itemHKL.Fc.getModule();
            sumFcFc += wHKL * Math.pow(itemHKL.Fc.getModule(), 2);
        }
        this.K = sumFoFc / sumFcFc;
    }

    public void improveK(double K) {
        this.K = K;
    }

    public void improvewExray(double wExray) {
        this.wExray = wExray;
    }

    public void improvewEcore(double wEcore) {
        this.wExray = wEcore;
    }

    public void calcEnergy(FragmentData FRAG) {
        /**
         * Calaculate Energy function value plus penalty function value.
         *
         * @param <code>Fragment Data</code>
         * @return Set corresponding instance values.
         */

        double partExray = 0;
        double partErest = 0;
        double sumFoFc = 0;
        double sumFcFc = 0;
        double sumFomkFc = 0;
        double sumwFomFcSq = 0;
        double sumwFoSq = 0;
        double sumFomFc = 0;
        double sumFo = 0;
        double sumFomFc_cond = 0;
        double sumFo_cond = 0;
        double Na = 0;
        double sumFcIm = 0;
        double sumFcRe = 0;
        double sumFomFcwk = 0;

        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            double wHKL = calcwHKL(itemHKL);
            sumFcIm = sumFcRe = 0;
            for (Fragment itemFrag : FRAG.getFragMass()) {
                ComplexNumber fScat = itemFrag.fragScattering(itemHKL, this.CELL, this.SYM);
                sumFcRe += fScat.getRe();
                sumFcIm += fScat.getIm();
            }
            itemHKL.Fc.setNum(sumFcRe, sumFcIm);
            sumFomkFc += Math.abs(itemHKL.Fo.getModule() - this.K * itemHKL.Fc.getModule());
            sumwFoSq += Math.abs(wHKL * Math.pow(itemHKL.Fo.getModule(), 2));
            sumwFomFcSq += Math.abs(wHKL * Math.pow(itemHKL.Fo.getModule() - itemHKL.Fc.getModule(), 2));
            sumFomFc_cond += (itemHKL.Fc.getModule() < itemHKL.Fo.getModule()) ? Math.abs(itemHKL.Fc.getModule() - itemHKL.Fo.getModule()) : 0.0;
            sumFo_cond += itemHKL.Fo.getModule();
            sumFomFc += Math.abs(itemHKL.Fc.getModule() - itemHKL.Fo.getModule());
            sumFo += itemHKL.Fo.getModule();
            sumFoFc += wHKL * itemHKL.Fo.getModule() * itemHKL.Fc.getModule();
            sumFcFc += wHKL * Math.pow(itemHKL.Fc.getModule(), 2);
            Na += wHKL * Math.pow(itemHKL.Fo.getModule(), 2);
            sumFomFcwk = wHKL * Math.pow(itemHKL.Fo.getModule() - this.K * itemHKL.Fc.getModule(), 2);
            itemHKL.Fc.setNum(sumFcRe, sumFcIm);
        }
        if (this.autoAdjustK) improveK(sumFoFc / sumFcFc);
        this.RI = sumFomkFc / sumFo;
        this.RII = sumFomFc_cond / sumFo_cond;
        this.RIII = Math.sqrt(sumwFomFcSq / sumwFoSq);
        this.RIV = sumFomFc / sumFo;

        if (this.OPT.contains("Xr1")) {
            partExray = sumFomFcwk;
        }
        if (this.OPT.contains("Xr2")) {
            partExray = sumFomFc;
        }
        if (this.OPT.contains("Xr3")) {
            partExray = sumFomkFc / Na;
        }

        if (this.OPT.contains("Xr4")) {
            partExray = sumwFomFcSq;
        }

        if (this.OPT.contains("Xr5")) {
            partExray = sumFomFc_cond;
        }

        if (this.OPT.contains("Pa")) {
            partExray = PattersonEnergy(this.SYM, FRAG);
        }

        this.Exray = partExray;

        if (this.OPT.contains("Re1")) {
            partErest += restrainFunction1(FRAG);
        }

        if (this.OPT.contains("Re2")) {
            partErest += restrainFunction2(FRAG);
        }

        this.Erest = partErest;
        this.Epenalty = this.PSI.Psi1(FRAG);
        this.Ecore = this.wExray * this.Exray + this.Erest;
        this.E = this.wEcore * this.Ecore + this.Epenalty;
    }
}

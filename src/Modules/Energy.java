package Modules;


import CrystTools.*;
import MathTools.FastMath;
import org.w3c.dom.DOMImplementation;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by Developer on 11.09.2015.
 */
public class Energy {
    private DiffractionData HKL;
    private UnitCell CELL;
    private Symmetry SYM;
    public String OPT;
    public double E;
    public double Exray;
    public double Erest;
    public double Ecore;
    public double Epenalty;
    public double wExray;
    public double wEcore;
    public double K;
    public double RI;
    public double RII;
    public double RIII;
    public double RIV;
    private PenaltyFunction PSI;
    private ArrayList<Double> statWhkl = new ArrayList<>();
    private ArrayList<Double> statFc = new ArrayList<>();


    public class aPair {
        public double Xi = 0;
        public double Yi = 0;
        public double Zi = 0;
        public double Xj = 0;
        public double Yj = 0;
        public double Zj = 0;
        public double DVdW = 0;
        public double R = 0;
        public aPair(double xi, double yi, double zi, double xj, double yj, double zj, double R, double DVdW) {
            Xi = xi;
            Yi = yi;
            Zi = zi;
            Xj = xj;
            Yj = yj;
            Zj = zj;
            this.R = R;
            this.DVdW = DVdW;
        }
    }


    public Energy(DiffractionData HKL, UnitCell CELL, Symmetry SYM, PenaltyFunction PSI) {
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
        this.K = 2.6891;
        this.RI = 0;
        this.RII = 0;
        this.RIII = 0;
        this.RIV = 0;
        calcwHKL(0);
    }

    public boolean ifInCell(double x, double y, double z){
        boolean condition = (x <= 1) && (x >= 0) && (y <= 1) && (y >= 0) && (z <= 1) && (z >= 0);
        if (condition) return true;
        return false;
    }

    public double restrainFunction(FragmentData FRAG){
        List<aPair> aPairs = genAtomsPairs(this.SYM, FRAG);
        double restrain = 0;
//        double c = 1;
//        double n = 2;
//        double m = 2;
//              if ((R != 0) && (R < itemApair.RVdW))
//                  partErest += Math.pow(Math.pow(c * RVdW, n) - Math.pow(R, n), m);
        for (aPair itemApair : aPairs) {
            if((itemApair.R != 0) && (itemApair.R < itemApair.DVdW)) {
                restrain += Math.pow(Math.min(itemApair.R - itemApair.DVdW, 0), 2);
            }
        }
            return restrain;
    }



    public List<aPair> genAtomsPairs(Symmetry SYM, FragmentData FRAG) {
        List<aPair> result = new ArrayList<aPair>();
        for (Fragment imFragment_i : FRAG.getFragMass()) {
            for (Fragment imFragment_j : FRAG.getFragMass()) {
                for (Atom imAtom_i : imFragment_i.getFragAtoms()) {
                    for (Atom imAtom_j : imFragment_j.getFragAtoms()) {
                        double DVdW = (imAtom_i.getAtomRVdW() + imAtom_j.getAtomRVdW());
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
                                double[] VECT_it = VECT_i;
                                for (int tX = -1; tX < 2; tX++){
                                    for (int tY = -1; tY < 2; tY++){
                                        for (int tZ = -1; tZ < 2; tZ++){
                                            double[] VECT_jt = FastMath.VpV(VECT_j, new double[]{tX, tY, tZ});
                                            if (ifInCell(
                                                    (VECT_jt[0] + VECT_it[0]) / 2.0,
                                                    (VECT_jt[1] + VECT_it[1]) / 2.0,
                                                    (VECT_jt[2] + VECT_it[2]) / 2.0)) {
                                                double R = this.CELL.calcDistance(VECT_it[0] - VECT_jt[0], VECT_it[1] - VECT_jt[1], VECT_it[2] - VECT_jt[2]);
                                                result.add(new aPair(
                                                        VECT_it[0], VECT_it[1], VECT_it[2],
                                                        VECT_jt[0], VECT_jt[1], VECT_jt[2], R, DVdW));
                                            }
                                        }
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
    public void calcwHKL(int mode){
        double FoMin = this.HKL.getHKL().get(0).Fsq;
        double FoMax = 0;
        double wHKL = 0;
        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            FoMin = Math.min(FoMin, itemHKL.Fsq);
            FoMax = Math.max(FoMax, itemHKL.Fsq);
        }
        double A = 2.0 * FoMin;
        double C = 2.0 / FoMax;
        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            switch (mode) {
                case 0:
                    wHKL = 1 / (A + itemHKL.Fsq + C * Math.pow(itemHKL.Fsq, 2));
                    break;
                case 1:
                    wHKL = 1 / Math.pow(itemHKL.Fsq, 2);
                    break;
                case 2:
                    wHKL = 1;
                    break;
                default:
                    break;
            }
            this.statWhkl.add(wHKL);
        }
    }


    public double PattersonFunction(double[] UVW){
        double V0 = this.CELL.getV();
        double sumPatt = 0;
        for (ReciprocalItem itemHKL : this.HKL.getHKL()) {
            sumPatt += itemHKL.Fsq * Math.cos(2 * Math.PI * ( itemHKL.h * UVW[0] + itemHKL.k * UVW[1] + itemHKL.l * UVW[2] ) );
        }
        return  sumPatt / V0;
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
        return  patt_E;
    }



    public void calcEnergy(FragmentData FRAG){
        /**
         * Calaculate Energy function value plus penalty function value.
         *
         * @param <code>Fragment Data</code>
         * @return Set corresponding instance values.
         */

        double partExray = 0;
        double partErest = 0;

        double[] cellTranslations = {this.CELL.getA(), this.CELL.getB(), this.CELL.getC()};
        Arrays.sort(cellTranslations);
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
        double Fc = 0;
        double sumFcIm = 0;
        double sumFcRe = 0;
        double Fo = 0;
        double sumFomFcwK = 0;

        for (int i = 0; i < this.HKL.getHKL().size(); i++) {
            sumFcIm = sumFcRe = 0;
            for (Fragment itemFrag : FRAG.getFragMass()) {
                double[] fScat = itemFrag.fragScattering(this.HKL.getHKL().get(i), this.CELL, this.SYM);
                sumFcRe +=  fScat[0];
                sumFcIm +=  fScat[1];
            }
            Fc = Math.sqrt(Math.pow(sumFcIm, 2) + Math.pow(sumFcRe, 2));
            Fo = Math.sqrt(this.HKL.getHKL().get(i).Fsq);
            sumFomkFc += Math.abs(Fo - this.K * Fc);
            sumwFoSq += Math.abs(this.statWhkl.get(i) * Math.pow(Fo, 2));
            sumwFomFcSq += Math.abs(this.statWhkl.get(i) * Math.pow(Fo - Fc, 2));
            if (Fc < Fo)  {
                sumFomFc_cond += Math.abs(Fc - Fo);
                sumFo_cond += Fo;
            }
            sumFomFc += Math.abs(Fc - Fo);
            sumFo += Fo;
            sumFoFc += this.statWhkl.get(i) * Fo * Fc;
            sumFcFc += this.statWhkl.get(i) * Math.pow(Fc, 2);
            Na += this.statWhkl.get(i) * Math.abs(this.HKL.getHKL().get(i).Fsq);
            sumFomFcwK = this.statWhkl.get(i) * Math.pow(Fo - this.K * Fc, 2);

        }
        this.K = sumFoFc / sumFcFc;
        this.RI = sumFomkFc / sumFo;
        this.RII = sumFomFc_cond / sumFo_cond;
        this.RIII = Math.sqrt(sumwFomFcSq / sumwFoSq);
        this.RIV = sumFomFc / sumFo;

        if (this.OPT.contains("Xr1")) {
            partExray = sumFomFcwK;
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

        if (this.OPT.contains("Pa")) {
            partExray = PattersonEnergy(this.SYM, FRAG);
        }

        this.Exray = partExray;

        if (this.OPT.contains("Re")) {
            partErest += restrainFunction(FRAG);
        }
        this.Erest = partErest;
        this.Epenalty = this.PSI.Psi1(FRAG);
        this.Ecore = this.wExray * this.Exray + this.Erest;
        this.E = this.wEcore * this.Ecore + this.Epenalty;
    }
}

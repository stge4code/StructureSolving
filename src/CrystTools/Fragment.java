package CrystTools;

import MathTools.FastMath;
import MathTools.SpecialFunction;

import java.io.*;
import java.util.ArrayList;

import static Utilities.ObjectsUtilities.deepClone;

/**
 * Created by Developer on 10.09.2015.
 */
public class Fragment implements Serializable {
    private int fragNum = 0;
    private double fragO = 0;
    private double fragU = 0;
    private  double fragMass = 0;
    private String fragScatType = "";
    private String fragOPT = "";

    private String fragName = "";

    private Double fragX = Double.valueOf(0);
    private Double fragY = Double.valueOf(0);
    private Double fragZ = Double.valueOf(0);

    private Double fragPhi1 = Double.valueOf(0);
    private Double fragPhi2 = Double.valueOf(0);
    private Double fragTheta = Double.valueOf(0);

    private ArrayList<Atom> fragAtoms = new ArrayList<>();
    private ArrayList<Atom> fragAtomsSource = new ArrayList<>();

    public String getFragName() {
        return fragName;
    }

    public Double getFragX() {
        return fragX;
    }

    public Double getFragY() {
        return fragY;
    }

    public Double getFragZ() {
        return fragZ;
    }

    public void setFragX(double fragX) {
        this.fragX = fragX;
    }

    public void setFragY(double fragY) {
        this.fragY = fragY;
    }

    public void setFragZ(double fragZ) {
        this.fragZ = fragZ;
    }

    public double getFragU() {
        return fragU;
    }

    public void setFragU(double fragU) {
        this.fragU = fragU;
    }

    public void setFragAtoms(ArrayList<Atom> fragAtoms) {
        this.fragAtoms = fragAtoms;
    }

    public Double getFragPhi1() {
        return fragPhi1;
    }

    public void setFragPhi1(double fragPhi1) {
        this.fragPhi1 = fragPhi1;
    }

    public Double getFragPhi2() {
        return fragPhi2;
    }

    public void setFragPhi2(double fragPhi2) {
        this.fragPhi2 = fragPhi2;
    }

    public Double getFragTheta() {
        return fragTheta;
    }

    public void setFragTheta(double fragTheta) {
        this.fragTheta = fragTheta;
    }

    public ArrayList<Atom> getFragAtoms() {
        return fragAtoms;
    }


    public void setFragO(double fragO) {
        this.fragO = fragO;
    }


    public int getFragNum() {
        return fragNum;
    }

    public void setFragNum(int fragNum) {
        this.fragNum = fragNum;
    }

    public double getFragO() {
        return fragO;
    }

    public String getFragOPT() {
        return fragOPT;
    }

    public void setFragOPT(String fragOPT) {
        this.fragOPT = fragOPT;
    }

    public String getFragScatType() {
        return fragScatType;
    }


    public void setFragScatType(String fragScatType) {
        this.fragScatType = fragScatType;
    }


    public Fragment(String fragName,
                    ArrayList<Atom> fragAtoms,
                    int fragNum,
                    double fragO,
                    double fragU,
                    String fragScatType,
                    String fragOPT) {
        this.fragName = fragName;
        this.fragAtoms = fragAtoms;
        this.fragNum = fragNum;
        this.fragO = fragO;
        this.fragU = fragU;
        this.fragScatType = fragScatType;
        this.fragOPT = fragOPT;
        fragFindMass();
        fragFindCenter(1);
        fragGenerateAtomsSource();
    }

    public void fragFindMass() {
        double fragMass = 0;
        for (Atom itemAtom : this.fragAtoms) {
            fragMass += itemAtom.getAtomM();
        }
        this.fragMass = fragMass;
    }

    public void fragFindCenter(int mode) {
        if (mode == 1) {
            double[] VECT = {0, 0, 0};
            for (Atom itemAtom : this.fragAtoms) {
                VECT[0] += itemAtom.getAtomX() * itemAtom.getAtomM() / this.fragMass;
                VECT[1] += itemAtom.getAtomY() * itemAtom.getAtomM() / this.fragMass;
                VECT[2] += itemAtom.getAtomZ() * itemAtom.getAtomM() / this.fragMass;
            }
            this.fragX = VECT[0];
            this.fragY = VECT[1];
            this.fragZ = VECT[2];
        } else if (mode == 2) {
            double[] VECT = {0, 0, 0};
            double atomsNum = this.fragAtoms.size();
            for(Atom itemAtom : this.fragAtoms){
                VECT[0] += itemAtom.getAtomX() / atomsNum;
                VECT[1] += itemAtom.getAtomY() / atomsNum;
                VECT[2] += itemAtom.getAtomZ() / atomsNum;
            }
            this.fragX = VECT[0];
            this.fragY = VECT[1];
            this.fragZ = VECT[2];
        }
    }


    public void fragGenerateAtomsSource() {

        this.fragAtomsSource = (ArrayList<Atom>) deepClone(fragAtoms);
        double[] VECT = {0, 0, 0};
        for (Atom itemAtom : this.fragAtomsSource) {
            VECT[0] += itemAtom.getAtomX() * itemAtom.getAtomM() / this.fragMass;
            VECT[1] += itemAtom.getAtomY() * itemAtom.getAtomM() / this.fragMass;
            VECT[2] += itemAtom.getAtomZ() * itemAtom.getAtomM() / this.fragMass;
        }
        for (Atom itemAtom : this.fragAtomsSource) {
            itemAtom.setAtomX(itemAtom.getAtomX() - VECT[0]);
            itemAtom.setAtomY(itemAtom.getAtomY() - VECT[1]);
            itemAtom.setAtomZ(itemAtom.getAtomZ() - VECT[2]);
        }
    }




    public void fragTranslate() {
        for (Atom itemAtomSource : this.fragAtomsSource) {
            double[] VECT = {itemAtomSource.getAtomX(), itemAtomSource.getAtomY(), itemAtomSource.getAtomZ()};
            VECT[0] += this.fragX;
            VECT[1] += this.fragY;
            VECT[2] += this.fragZ;
            Atom itemAtom = this.fragAtoms.get(this.fragAtomsSource.indexOf(itemAtomSource));
            itemAtom.setAtomX(VECT[0]);
            itemAtom.setAtomY(VECT[1]);
            itemAtom.setAtomZ(VECT[2]);
        }

    }

    public void fragModifyParameters(UnitCell CELL,
                                     double dX, double dY, double dZ,
                                     double dPhi1, double dPhi2, double dTheta,
                                     double dO, double dU) {

        this.fragO += dO;
        this.fragU += dU;
        this.fragX += dX;
        this.fragY += dY;
        this.fragZ += dZ;
        this.fragPhi1 += dPhi1;
        this.fragPhi2 += dPhi2;
        this.fragTheta += dTheta;

        fragTranslate();
        fragRotate(CELL);
    }


    public void fragRotate(UnitCell CELL) {
        double scaleFactor =  2 * Math.PI;

        double[][] M = new double[3][3];
        double cT = Math.cos(scaleFactor * this.fragTheta);
        double sT = Math.sin(scaleFactor * this.fragTheta);
        double cP1 = Math.cos(scaleFactor * this.fragPhi1);
        double sP1 = Math.sin(scaleFactor * this.fragPhi1);
        double cP2 = Math.cos(scaleFactor * this.fragPhi2);
        double sP2 = Math.sin(scaleFactor * this.fragPhi2);



        M[0][0] = cP1 * cP2 - cT * sP1 * sP2;
        M[0][1] = -cP1 * sP2 - cT * sP1 * cP2;
        M[0][2] = sP1 * sT;

        M[1][0] = sP1 * cP2 + cT * cP1 * sP2;
        M[1][1] = -sP1 * sP2 + cT * cP1 * cP2;
        M[1][2] = -cP1 * sT;

        M[2][0] = sP2 * sT;
        M[2][1] = cP2 * sT;
        M[2][2] = cT;



        for (Atom itemAtomSource : this.fragAtomsSource) {
            double[] VECT = {itemAtomSource.getAtomX(), itemAtomSource.getAtomY(), itemAtomSource.getAtomZ()};
            VECT = CELL.f2c(VECT);
            VECT = FastMath.MmV(M, VECT);
            VECT = CELL.c2f(VECT);
            VECT[0] += this.fragX;
            VECT[1] += this.fragY;
            VECT[2] += this.fragZ;
            Atom itemAtom = this.fragAtoms.get(this.fragAtomsSource.indexOf(itemAtomSource));
            itemAtom.setAtomX(VECT[0]);
            itemAtom.setAtomY(VECT[1]);
            itemAtom.setAtomZ(VECT[2]);
        }
    }


    public void fragParametersToOrder() {
        this.fragO = (this.fragO < 0) ? 1.0 + this.fragO % 1.0 : ((this.fragO % 1.0 == 0) ? 1.0 : this.fragO % 1.0);
        this.fragX = (this.fragX < 0) ? 1.0 + this.fragX % 1.0 : this.fragX % 1.0;
        this.fragY = (this.fragY < 0) ? 1.0 + this.fragY % 1.0 : this.fragY % 1.0;
        this.fragZ = (this.fragZ < 0) ? 1.0 + this.fragZ % 1.0 : this.fragZ % 1.0;
        //this.fragPhi1 = (this.fragPhi1 < 0) ? 1.0 + this.fragPhi1 % 1.0 : this.fragPhi1 % 1.0;
        //this.fragPhi2 = (this.fragPhi2 < 0) ? 1.0 + this.fragPhi2 % 1.0 : this.fragPhi2 % 1.0;
        //this.fragTheta = (this.fragTheta < 0) ? 1.0 + this.fragTheta % 1.0 : this.fragTheta % 1.0;
    }




    public void fragApplySymmetry(double[][] fragSYMr, double[] fragSYMt) {
        for (Atom itemAtom : this.fragAtoms) {
            double[] VECT = {itemAtom.getAtomX(), itemAtom.getAtomY(), itemAtom.getAtomZ()};
            VECT = FastMath.VpV(fragSYMt, FastMath.MmV(fragSYMr, VECT));
            itemAtom.setAtomX(VECT[0]);
            itemAtom.setAtomY(VECT[1]);
            itemAtom.setAtomZ(VECT[2]);
        }
    }


    public double[] fragScattering(ReciprocalItem hkl, UnitCell cell, Symmetry sym) {
        return fragScattering(hkl, cell, sym, this.fragO, this.fragScatType);
    }

    public double[] fragScattering(ReciprocalItem HKL,
                                   UnitCell CELL,
                                   Symmetry SYM,
                                   double fragOccupancy,
                                   String scatType) {
        double[] Fhkl = {0.0, 0.0};
        double atomScattering = 0;
        double argPhase = 0;
        double ScatVect = HKL.scatvect;
        double ScatVectSq = Math.pow(HKL.scatvect, 2);
        double Debye_Coefficient = Math.exp(-8.0 * Math.pow(Math.PI, 2) * this.fragU * ScatVectSq);
        if (scatType.contains("F1")) {
            atomScattering = fragOccupancy * Debye_Coefficient * this.fragAtoms.get(0).atomScattering(ScatVectSq);
            double D = fragDiameter(CELL);
            double[] VECTCORD = {this.fragX, this.fragY, this.fragZ};
            for (SymmetryItem itemSYM : SYM.getSymMass()) {
                double[] VECT = FastMath.VpV(itemSYM.getSymT(), FastMath.MmV(itemSYM.getSymR(), VECTCORD));
                argPhase = 2 * Math.PI * (VECT[0] * HKL.h + VECT[1] * HKL.k + VECT[2] * HKL.l);
                Fhkl[0] += Math.cos(argPhase);
                if (SYM.getLATT() < 0) {
                    Fhkl[1] += Math.sin(argPhase);
                } else {
                    Fhkl[1] = 0;
                }
            }
            double BessScat = SpecialFunction.j0(2 * Math.PI * ScatVect * D);
            Fhkl[0] *= 60 * atomScattering * BessScat;
            Fhkl[1] *= 60 * atomScattering * BessScat;
        } else if (scatType.contains("F0")) {
            for (Atom itemAtom : this.fragAtoms) {
                atomScattering = itemAtom.atomScattering(ScatVectSq);
                double[] VECTCORD = {itemAtom.getAtomX(), itemAtom.getAtomY(), itemAtom.getAtomZ()};
                for (SymmetryItem itemSYM : SYM.getSymMass()) {
                    double[] VECT = FastMath.VpV(itemSYM.getSymT(), FastMath.MmV(itemSYM.getSymR(), VECTCORD));
                    argPhase = 2 * Math.PI * (VECT[0] * HKL.h + VECT[1] * HKL.k + VECT[2] * HKL.l);
                    Fhkl[0] += atomScattering * Math.cos(argPhase);
                    Fhkl[1] +=  atomScattering * Math.sin(argPhase);
                    if (SYM.getLATT() < 0) {
                        Fhkl[1] +=  atomScattering * Math.sin(argPhase);
                    } else {
                        Fhkl[1] = 0;
                    }
                }

            }
            Fhkl[0] *= fragOccupancy * Debye_Coefficient;
            Fhkl[1] *= fragOccupancy * Debye_Coefficient;
        }
        return Fhkl;
    }

    public double fragDiameter(UnitCell CELL) {
        double D = 0;
        int numAtoms = this.fragAtoms.size();
        ArrayList<Double> rMass = new ArrayList();
        for (Atom itemAtom : this.fragAtoms) {
            rMass.add(CELL.calcDistance(
                    itemAtom.getAtomX() - this.fragX,
                    itemAtom.getAtomY() - this.fragY,
                    itemAtom.getAtomZ() - this.fragZ));
        }
        for (Double itemR : rMass) {
            D += 2.0 / numAtoms * itemR;
        }
        return D;
    }

}
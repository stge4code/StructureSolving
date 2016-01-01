package CrystTools;

/**
 * Created by Developer on 10.09.2015.
 */
public class Atom implements java.io.Serializable {
//    public Atom(String atomName, double atomX, double atomY, double atomZ, double[] atomScat, double atomRVdW, double atomU) {
//        super(atomName, atomX, atomY, atomZ, atomScat, atomRVdW, atomU);
//    }

    private String atomName = "";
    private double atomX = 0;
    private double atomY = 0;
    private double atomZ = 0;
    private double[] atomScat = new double[9];
    private double atomRVdW = 0;
    private double atomM = 0;


    public String getAtomName() {
        return atomName;
    }

    public void setAtomName(String atomName) {
        this.atomName = atomName;
    }

    public double getAtomX() {
        return atomX;
    }

    public void setAtomX(double atomX) {
        this.atomX = atomX;
    }

    public double getAtomY() {
        return atomY;
    }

    public void setAtomY(double atomY) {
        this.atomY = atomY;
    }

    public double getAtomZ() {
        return atomZ;
    }

    public void setAtomZ(double atomZ) {
        this.atomZ = atomZ;
    }

    public double[] getAtomScat() {
        return atomScat;
    }

    public void setAtomScat(double[] atomScat) {
        this.atomScat = atomScat;
    }

    public double getAtomRVdW() {
        return atomRVdW;
    }

    public double getAtomM() {
        return atomM;
    }

    public void setAtomM(double atomM) {
        this.atomM = atomM;
    }

    public void setAtomRVdW(double atomRVdW) {
        this.atomRVdW = atomRVdW;
    }


    public Atom(String atomName,
                double atomX,
                double atomY,
                double atomZ,
                double[] atomScat,
                double atomRVdW,
                double atomM) {
        this.atomName = atomName;
        this.atomX = atomX;
        this.atomY = atomY;
        this.atomZ = atomZ;
        this.atomScat = atomScat;
        this.atomRVdW = atomRVdW;
        this.atomM = atomM;

    }

    public double atomScattering(double scatVectSq) {
        //double scatVectSq = Math.pow(cell.calcScatVect(hkl), 2);
        double sumExps = this.atomScat[8];
        for (int i = 0; i < 4; i += 2) {
            sumExps += this.atomScat[i] * Math.exp(-this.atomScat[i + 1] * scatVectSq);
        }
        return sumExps;
    }
}

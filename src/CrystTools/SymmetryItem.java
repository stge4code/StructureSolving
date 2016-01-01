package CrystTools;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Created by Administrator on 14.09.2015.
 */
public class SymmetryItem implements Serializable {
    private double[][] symR = new double[3][3];
    private double[] symT = new double[3];

    public double[][] getSymR() {
        return symR;
    }

    public void setSymR(double[][] symM) {
        this.symR = symM;
    }

    public double[] getSymT() {
        return symT;
    }

    public void setSymT(double[] symT) {
        this.symT = symT;
    }

    public SymmetryItem(double[][] symR, double[] symT) {
        this.symR = symR;
        this.symT = symT;
    }
}

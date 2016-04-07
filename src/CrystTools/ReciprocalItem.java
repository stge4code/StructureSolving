package CrystTools;

import MathTools.ComplexNumber;
import MathTools.FastMath;

/**
 * Created by Developer on 10.09.2015.
 */
public class ReciprocalItem {
    public double h;
    public double k;
    public double l;
    public double I;
    public double sigmaI;
    public double batchNumber;
    public double scatvect;
    public ComplexNumber Fo = new ComplexNumber();
    public ComplexNumber Fc = new ComplexNumber();

    public ReciprocalItem(double h, double k, double l, double i, double sigmai, double scatvect, double batchnumber) {
        this.h = h;
        this.k = k;
        this.l = l;
        this.I = i;
        this.sigmaI = sigmai;
        this.scatvect = scatvect;
        this.batchNumber = batchnumber;
        updateFoModule(i);
    }

    public ReciprocalItem(ReciprocalItem other) {
        this.h = other.h;
        this.k = other.k;
        this.l = other.l;
        this.I = other.I;
        this.sigmaI = other.sigmaI;
        this.scatvect = other.scatvect;
        this.batchNumber = other.batchNumber;
        this.Fc.setNum(other.Fc);
        this.Fo.setNum(other.Fo);
        updateFoModule(other.I);
    }

    public ReciprocalItem modify(String hcond, String kcond, String lcond) {
        ReciprocalItem result = new ReciprocalItem(this);
        result.h = (int) FastMath.eval(hcond
                .replaceAll("h", "(" + Double.toString(this.h) + ")")
                .replaceAll("k", "(" + Double.toString(this.k) + ")")
                .replaceAll("l", "(" + Double.toString(this.l) + ")"));
        result.k = (int) FastMath.eval(kcond
                .replaceAll("h", "(" + Double.toString(this.h) + ")")
                .replaceAll("k", "(" + Double.toString(this.k) + ")")
                .replaceAll("l", "(" + Double.toString(this.l) + ")"));
        result.l = (int) FastMath.eval(lcond
                .replaceAll("h", "(" + Double.toString(this.h) + ")")
                .replaceAll("k", "(" + Double.toString(this.k) + ")")
                .replaceAll("l", "(" + Double.toString(this.l) + ")"));
        return result;
    }

    public ReciprocalItem modify(double hmodified, double kmodified, double lmodified) {
        ReciprocalItem result = new ReciprocalItem(this);
        result.h = hmodified;
        result.k = kmodified;
        result.l = lmodified;
        return result;
    }

    public static boolean compare(ReciprocalItem itemi, ReciprocalItem itemj) {
        return ((itemi.h == itemj.h) && (itemi.k == itemj.k) && (itemi.l == itemj.l));
    }




    private void updateFoModule(double i) {
        if (i >= 0) this.Fo.setModule(Math.sqrt(i));
    }
}

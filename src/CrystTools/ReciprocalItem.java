package CrystTools;

import MathTools.ComplexNumber;
import MathTools.FastMath;

/**
 * Created by Developer on 10.09.2015.
 */
public class ReciprocalItem {
    public int h;
    public int k;
    public int l;
    public double qh;
    public double qk;
    public double ql;
    public double I;
    public double sigmaI;
    public double batchNumber;
    public double scatvect;
    public ComplexNumber Fo = new ComplexNumber();
    public ComplexNumber Fc = new ComplexNumber();

    public ReciprocalItem(int h, int k, int l, double i, double sigmai, double scatvect, double batchnumber, double qh, double qk, double ql) {
        this.h = h;
        this.k = k;
        this.l = l;
        this.I = i;
        this.sigmaI = sigmai;
        this.scatvect = scatvect;
        this.batchNumber = batchnumber;
        updateFoModule(i);
        this.qh = qh;
        this.qk = qk;
        this.ql = ql;
    }

    public ReciprocalItem(ReciprocalItem other) {
        this.h = other.h;
        this.k = other.k;
        this.l = other.l;
        this.qh = other.qh;
        this.qk = other.qk;
        this.ql = other.ql;
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
                .replaceAll("h", "(" + Integer.toString(this.h) + ")")
                .replaceAll("k", "(" + Integer.toString(this.k) + ")")
                .replaceAll("l", "(" + Integer.toString(this.l) + ")"));
        result.k = (int) FastMath.eval(kcond
                .replaceAll("h", "(" + Integer.toString(this.h) + ")")
                .replaceAll("k", "(" + Integer.toString(this.k) + ")")
                .replaceAll("l", "(" + Integer.toString(this.l) + ")"));
        result.l = (int) FastMath.eval(lcond
                .replaceAll("h", "(" + Integer.toString(this.h) + ")")
                .replaceAll("k", "(" + Integer.toString(this.k) + ")")
                .replaceAll("l", "(" + Integer.toString(this.l) + ")"));
        return result;
    }

    public ReciprocalItem modify(int hmodified, int kmodified, int lmodified) {
        ReciprocalItem result = new ReciprocalItem(this);
        result.h = hmodified;
        result.k = kmodified;
        result.l = lmodified;
        return result;
    }

    public static boolean compare(ReciprocalItem itemi, ReciprocalItem itemj) {
        return ((itemi.h == itemj.h) &&
                (itemi.k == itemj.k) &&
                (itemi.l == itemj.l) &&
                (itemi.qh == itemj.qh) &&
                (itemi.qk == itemj.qk) &&
                (itemi.ql == itemj.ql));
    }




    private void updateFoModule(double i) {
        if (i >= 0) this.Fo.setModule(Math.sqrt(i));
    }
}

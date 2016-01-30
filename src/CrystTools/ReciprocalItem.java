package CrystTools;

import MathTools.ComplexNumber;

/**
 * Created by Developer on 10.09.2015.
 */
public class ReciprocalItem {
    public int h;
    public int k;
    public int l;
    public double I;
    public double sigmaI;
    public double batchNumber;
    public double scatvect;
    public ComplexNumber Fo = new ComplexNumber();
    public ComplexNumber Fc = new ComplexNumber();

    public ReciprocalItem(int h, int k, int l, double i, double sigmai, double scatvect, double batchnumber) {
        this.h = h;
        this.k = k;
        this.l = l;
        this.I = i;
        this.sigmaI = sigmai;
        this.scatvect = scatvect;
        this.batchNumber = batchnumber;
        this.Fo.setModule(Math.sqrt(i));
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
    }
}

package MathTools;

/**
 * Created by Developer on 28.01.2016.
 */
public class ComplexNumber {
    private double Im;
    private double Re;
    private double module;
    private double phase;

    public ComplexNumber(double re, double im) {
        this.Re = re;
        this.Im = im;
        this.module = Math.sqrt(Math.pow(re, 2) + Math.pow(im, 2));
        this.phase = (this.module != 0) ? Math.atan(im / re) : 0.0;
    }

    public ComplexNumber(double [] reim) {
        this.Re = reim[0];
        this.Im = reim[1];
        this.module = Math.pow(reim[0], 2) + Math.pow(reim[1], 2);
        this.phase = (this.module != 0) ? Math.atan(reim[1] / reim[0]) : 0.0;
    }

    public double addToIm(double im) {
        this.Im += im;
        return this.Im;
    }

    public double addToRe(double re) {
        this.Re += re;
        return this.Re;
    }

    public double getIm() {
        return Im;
    }

    public void setIm(double im) {
        Im = im;
    }

    public double getRe() {
        return Re;
    }

    public void setRe(double re) {
        Re = re;
    }

    public double getModule() {
        return module;
    }

    public void setModule(double module) {
        this.module = module;
    }

    public double getPhase() {
        return phase;
    }

    public void setPhase(double phase) {
        this.phase = phase;
    }
}

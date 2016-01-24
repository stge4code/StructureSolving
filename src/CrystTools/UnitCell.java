package CrystTools;

import MathTools.FastMath;
import Utilities.ObjectsUtilities;

import java.io.*;
import java.net.SocketPermission;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by Administrator on 11.09.2015.
 */
public class UnitCell implements Serializable {
    private String unitCellFileName;

    private double a;
    private double b;
    private double c;
    private double alpha;
    private double beta;
    private double gamma;
    private double V;

    private double aStar;
    private double bStar;
    private double cStar;
    private double alphaStar;
    private double betaStar;
    private double gammaStar;
    private double Vstar;

    private double[][] Gstar;
    private double[][] G;


    private double cA;
    private double sA;
    private double cB;
    private double sB;
    private double cG;
    private double sG;

    private double cAstar;
    private double sAstar;
    private double cBstar;
    private double sBstar;
    private double cGstar;
    private double sGstar;


    public void setA(double a) {
        this.a = a;
    }

    public void setB(double b) {
        this.b = b;
    }

    public void setC(double c) {
        this.c = c;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    public void setGamma(double gamma) {
        this.gamma = gamma;
    }

    public double getV() {
        return V;
    }

    public void setV(double v) {
        V = v;
    }

    public double getA() {
        return a;
    }

    public double getB() {
        return b;
    }

    public double getC() {
        return c;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getBeta() {
        return beta;
    }

    public double getGamma() {
        return gamma;
    }

    public void calcVolume() {


        this.V = this.a * this.b * this.c *
                Math.sqrt(1 - Math.pow(this.cA, 2.0) - Math.pow(this.cB, 2.0) - Math.pow(this.cG, 2.0) + 2 * this.cA * this.cB * this.cG);
    }


    public void calcReciprocalVect() {

        this.aStar = Math.abs(this.b * this.c * FastMath.d2rSin(this.alpha)) / this.V;
        this.bStar = Math.abs(this.a * this.c * FastMath.d2rSin(this.beta)) / this.V;
        this.cStar = Math.abs(this.a * this.b * FastMath.d2rSin(this.gamma)) / this.V;

//        this.alphaStar = Math.asin(this.V / (this.a * this.b * this.c *
//                FastMath.d2rCos(this.beta) * FastMath.d2rCos(this.gamma))) * 180 / Math.PI;


        this.alphaStar = Math.acos((this.cB * this.cG - this.cA) / (this.sB * this.sG)) * 180 / Math.PI;
        this.betaStar = Math.acos((this.cA * this.cG - this.cB) / (this.sA * this.sG)) * 180 / Math.PI;
        this.gammaStar = Math.acos((this.cA * this.cB - this.cG) / (this.sA * this.sB)) * 180 / Math.PI;

//        this.betaStar = Math.asin(this.V / (this.a * this.b * this.c *
//                FastMath.d2rCos(this.alpha) * FastMath.d2rCos(this.gamma))) * 180 / Math.PI;
//        this.gammaStar = Math.asin(this.V / (this.a * this.b * this.c *
//                FastMath.d2rCos(this.alpha) * FastMath.d2rCos(this.beta))) * 180 / Math.PI;

        this.Vstar = 1 / this.V;


    }

    public void calcG() {
        double[][] G = new double[3][3];
        G[0][0] = this.a * this.a;
        G[1][1] = this.b * this.b;
        G[2][2] = this.c * this.c;

        G[0][1] = this.a * this.b * cG;
        G[0][2] = this.a * this.c * cB;

        G[1][0] = this.a * this.b * cG;
        G[1][2] = this.b * this.c * cA;

        G[2][0] = this.a * this.c * cB;
        G[2][1] = this.b * this.c * cA;
        this.G = G;
    }

    public void calcGstar() {
        double[][] Gstar = new double[3][3];

        Gstar[0][0] = this.aStar * this.aStar;
        Gstar[1][1] = this.bStar * this.bStar;
        Gstar[2][2] = this.cStar * this.cStar;

        Gstar[0][1] = this.aStar * this.bStar * cGstar;
        Gstar[0][2] = this.aStar * this.cStar * cBstar;

        Gstar[1][0] = this.aStar * this.bStar * cGstar;
        Gstar[1][2] = this.bStar * this.cStar * cAstar;

        Gstar[2][0] = this.aStar * this.cStar * cBstar;
        Gstar[2][1] = this.bStar * this.cStar * cAstar;

        this.Gstar = Gstar;
    }

    public double calcDistance(double[] VECT) {
        return calcDistance(VECT[0], VECT[1], VECT[2]);
    }

    public double calcDistance(double x, double y, double z) {
        double[] VECT = {x, y, z};
        return Math.sqrt(FastMath.VmV(VECT, FastMath.MmV(this.G, VECT)));
    }

    public double[] c2f(double[] VECT) {
        return c2f(VECT[0], VECT[1], VECT[2]);
    }

    public double[] c2f(double x, double y, double z) {
        double[] VECT = new double[3];

        //vect[0] = 1 / this.a * x - cG / this.a / sG * y + (cA * cG - cB) / (sG * this.V / this.b / this.c) * z;
        //vect[1] = 1 / (this.b * sG) * y + (cB * cG - cA) / (sG * this.V / this.a / this.c) * z;
        //vect[2] = sG / (this.V / this.b / this.a) * z;


        VECT[0] = 1 / this.a * x - this.cG / this.a / this.sG * y +
                (this.b * this.c * this.cG * (this.cA - this.cB * this.cG) / this.sG - this.b * this.c * this.cB * this.sG) / this.V * z;
        VECT[1] = 1 / this.b / this.sG * y - this.a * this.c * (this.cA - this.cB * this.cG) / this.sG / this.V * z;
        VECT[2] = this.a * this.b * this.sG / this.V * z;

        return VECT;
    }

    public double[] f2c(double[] VECT) {
        return f2c(VECT[0], VECT[1], VECT[2]);
    }

    public double[] f2c(double x, double y, double z) {

        double[] VECT = new double[3];

        VECT[0] = this.a * x + this.b * this.cG * y + this.c * this.cB * z;
        VECT[1] = this.b * this.sG * y + this.c * (this.cA - this.cB * this.cG) / this.sG * z;
        VECT[2] = this.V / this.sG / this.a / this.b * z;

        return VECT;
    }


    /**
     * This methos is used for calculation of Scattering vector: S = 1 / (2 * d) = sin(theta) / lambda;
     *
     * @param hkl
     * @return
     */
    public double calcScatVect(ReciprocalItem hkl) {
        return calcScatVect(hkl.h, hkl.k, hkl.l);
    }

    public double calcScatVect(int h, int k, int l) {
        /*
        double h = hkl.h;
        double k = hkl.k;
        double l = hkl.l;
        double[][] S = new double[3][3];
        S[0][0] = Math.pow(this.b * this.c * FastMath.d2rSin(this.alpha), 2.0);
        S[1][1] = Math.pow(this.a * this.c * FastMath.d2rSin(this.beta), 2.0);
        S[2][2] = Math.pow(this.a * this.b * FastMath.d2rSin(this.gamma), 2.0);
        S[0][1] =  this.a * this.b * Math.pow(this.c, 2.0) *
                (FastMath.d2rCos(this.alpha) * FastMath.d2rCos(this.beta) - FastMath.d2rCos(this.gamma));
        S[1][2] =  this.c * this.b * Math.pow(this.a, 2.0) *
                (FastMath.d2rCos(this.beta) * FastMath.d2rCos(this.gamma) - FastMath.d2rCos(this.alpha));
        S[0][2] =  this.a * this.c * Math.pow(this.b, 2.0) *
                (FastMath.d2rCos(this.gamma) * FastMath.d2rCos(this.alpha) - FastMath.d2rCos(this.beta));
        return 2 * .5 * Math.sqrt(S[0][0] * h * h +
                S[1][1] * k * k +
                S[2][2] * l * l +
                2 * (S[0][1] * h * k + S[1][2] * k * l + S[0][2] * h * l)) / this.V;
        */

        double[] Vhkl = {h, k, l};
        return Math.sqrt(FastMath.VmV(Vhkl, FastMath.MmV(this.Gstar, Vhkl))) / 2.0;
    }


    public void printCell() {
        System.out.printf("\nCELL = % 10.3f % 10.3f % 10.3f\n       % 10.3f % 10.3f % 10.3f\n       % 10.3f\n",
                this.a, this.b, this.c,
                this.alpha, this.beta, this.gamma,
                this.V);
    }

    public void printCellstar() {
        System.out.printf("\nCELL* = % 10.3f % 10.3f % 10.3f\n        % 10.3f % 10.3f % 10.3f\n        % 10.3f\n",
                this.aStar, this.bStar, this.cStar,
                this.alphaStar, this.betaStar, this.gammaStar,
                this.Vstar);
    }

    public void printG() {
        System.out.printf("\n    % 10.3f % 10.3f % 10.3f\nG = % 10.3f % 10.3f % 10.3f\n    % 10.3f % 10.3f % 10.3f\n",
                G[0][0], G[0][1], G[0][2],
                G[1][0], G[1][1], G[1][2],
                G[2][0], G[2][1], G[2][2]);
    }

    public void printGstar() {
        System.out.printf("\n     % 10.3f % 10.3f % 10.3f\nG* = % 10.3f % 10.3f % 10.3f\n     % 10.3f % 10.3f % 10.3f\n",
                Gstar[0][0], Gstar[0][1], Gstar[0][2],
                Gstar[1][0], Gstar[1][1], Gstar[1][2],
                Gstar[2][0], Gstar[2][1], Gstar[2][2]);
    }

    public void printGmGstar() {
        double[][] GmGstar = FastMath.MmM(G, Gstar);
        System.out.printf("\n      % 10.3f % 10.3f % 10.3f\nGG* = % 10.3f % 10.3f % 10.3f\n      % 10.3f % 10.3f % 10.3f\n",
                GmGstar[0][0], GmGstar[0][1], GmGstar[0][2],
                GmGstar[1][0], GmGstar[1][1], GmGstar[1][2],
                GmGstar[2][0], GmGstar[2][1], GmGstar[2][2]);
    }

    public void calcTrigonometry() {
        this.cA = FastMath.d2rCos(this.alpha);
        this.sA = FastMath.d2rSin(this.alpha);
        this.cB = FastMath.d2rCos(this.beta);
        this.sB = FastMath.d2rSin(this.beta);
        this.cG = FastMath.d2rCos(this.gamma);
        this.sG = FastMath.d2rSin(this.gamma);
    }

    public void calcTrigonometryStar() {
        this.cAstar = FastMath.d2rCos(this.alphaStar);
        this.sAstar = FastMath.d2rSin(this.alphaStar);
        this.cBstar = FastMath.d2rCos(this.betaStar);
        this.sBstar = FastMath.d2rSin(this.betaStar);
        this.cGstar = FastMath.d2rCos(this.gammaStar);
        this.sGstar = FastMath.d2rSin(this.gammaStar);

    }

    public UnitCell(String unitCellFileName) {
        this.unitCellFileName = unitCellFileName;
        List<String> input = ObjectsUtilities.getContentFromFile(this.unitCellFileName);
        for (String s : input) {
            if (s.indexOf("CELL") != -1) {
                try {
                    ArrayList<String> allMatches = new ArrayList<String>();
                    Matcher m = Pattern.compile("[-+]?[\\d]*\\.?[\\d]+").matcher(s.substring(4));
                    while (m.find()) {
                        allMatches.add(m.group());
                    }
                    this.a = (Double.valueOf(allMatches.get(1)).doubleValue());
                    this.b = (Double.valueOf(allMatches.get(2)).doubleValue());
                    this.c = (Double.valueOf(allMatches.get(3)).doubleValue());
                    this.alpha = (Double.valueOf(allMatches.get(4)).doubleValue());
                    this.beta = (Double.valueOf(allMatches.get(5)).doubleValue());
                    this.gamma = (Double.valueOf(allMatches.get(6)).doubleValue());
                    this.calcTrigonometry();
                    this.calcVolume();
                    this.calcReciprocalVect();
                    this.calcTrigonometryStar();
                    this.calcG();
                    this.calcGstar();
                } catch (NumberFormatException e) {
                }
            }

        }
    }
}

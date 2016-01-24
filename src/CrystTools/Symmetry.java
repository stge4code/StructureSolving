package CrystTools;

import MathTools.FastMath;
import Utilities.ObjectsUtilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static Utilities.ObjectsUtilities.deepClone;

/**
 * Created by Administrator on 14.09.2015.
 */
public class Symmetry {
    private List<SymmetryItem> symMass = new ArrayList<>();
    private String symFileName;
    private int LATT = 0;

    public List<SymmetryItem> getSymMass() {
        return symMass;
    }

    public void setSymMass(ArrayList<SymmetryItem> symMass) {
        this.symMass = symMass;
    }

    public Symmetry(String symFileName) {
        List<SymmetryItem> symMassSourse = new ArrayList<SymmetryItem>();
        this.symFileName = symFileName;
        List<String> input = ObjectsUtilities.getContentFromFile(this.symFileName);
        for (String s : input) {
                    if(s.indexOf("LATT")!= -1) {
                        this.LATT = Integer.parseInt(s.substring(4).replaceAll("\\s+", ""));
                    }

                    if(s.indexOf("SYMM")!= -1) {
                        try {

                            List<String> allMatches = new ArrayList<String>();
                            Matcher m = Pattern.compile("[^,]+").matcher(s.substring(4).replaceAll("\\s+", ""));
                            while (m.find()) {
                                allMatches.add(m.group());
                            }

                            double[][] SYMr = new double[3][3];
                            double[] SYMt = new double[3];
                            String sTemp = "";

                            for(int i = 0; i < 3; i++){
                                sTemp = allMatches.get(i);
                                SYMt[i] = FastMath.eval(sTemp
                                        .replaceAll("X", "0.0")
                                        .replaceAll("Y", "0.0")
                                        .replaceAll("Z", "0.0"));
                            }
                            for(int i = 0; i < 3; i++){
                                sTemp = allMatches.get(i);
                                SYMr[i][0] = -SYMt[i] + FastMath.eval(sTemp
                                        .replaceAll("X", "1.0")
                                        .replaceAll("Y", "0.0")
                                        .replaceAll("Z", "0.0"));
                                SYMr[i][1] = -SYMt[i] + FastMath.eval(sTemp
                                        .replaceAll("X", "0.0")
                                        .replaceAll("Y", "1.0")
                                        .replaceAll("Z", "0.0"));
                                SYMr[i][2] = -SYMt[i] + FastMath.eval(sTemp
                                        .replaceAll("X", "0.0")
                                        .replaceAll("Y", "0.0")
                                        .replaceAll("Z", "1.0"));
                            }
                            this.symMass.add(new SymmetryItem(SYMr, SYMt));
                            if (this.LATT > 0) {
                                SYMt = (double[]) deepClone(SYMt);
                                SYMr = (double[][]) deepClone(SYMr);
                                for(int i = 0; i < 3; i++){
                                    SYMt[i] *= -1;
                                   for (int j = 0; j < 3; j++) {
                                       SYMr[i][j] *= -1;
                                   }
                                }
                                this.symMass.add(new SymmetryItem(SYMr, SYMt));
                            }
                            SYMr = null;
                            SYMt = null;
                            allMatches.clear();
                        } catch (NumberFormatException s2nSYM) {
                        }


                    }
                }

                double[][] SYMr = new double[3][3];
                double[] SYMt = new double[3];
                SYMt[0] = SYMt[1] = SYMt[2] = 0;
                SYMr[0][0] = SYMr[1][1] = SYMr[2][2] = 1;
                SYMr[0][1] = SYMr[0][2] = SYMr[1][0] = SYMr[1][2] = SYMr[2][0] = SYMr[2][1] = 0;
                this.symMass.add(new SymmetryItem(SYMr, SYMt));
                if (LATT > 0) {
                    SYMt = (double[]) deepClone(SYMt);
                    SYMr = (double[][]) deepClone(SYMr);
                    for(int i = 0; i < 3; i++){
                        SYMt[i] *= -1;
                        for (int j = 0; j < 3; j++) {
                            SYMr[i][j] *= -1;
                        }
                    }
                    this.symMass.add(new SymmetryItem(SYMr, SYMt));
                }
                SYMr = null;
                SYMt = null;
                LATTirazation();


    }

    public int getLATT() {
        return LATT;
    }

    public void setLATT(int LATT) {
        this.LATT = LATT;
    }

    public void printSymmetry(){
        System.out.printf("\nSymmetry operations:\n");
        for(SymmetryItem itemSym : this.symMass){
            System.out.printf("\n% 7.3f % 7.3f % 7.3f    % 7.3f\n% 7.3f % 7.3f % 7.3f  + % 7.3f\n% 7.3f % 7.3f % 7.3f    % 7.3f\n",
                    itemSym.getSymR()[0][0],
                    itemSym.getSymR()[0][1],
                    itemSym.getSymR()[0][2],
                    itemSym.getSymT()[0],
                    itemSym.getSymR()[1][0],
                    itemSym.getSymR()[1][1],
                    itemSym.getSymR()[1][2],
                    itemSym.getSymT()[1],
                    itemSym.getSymR()[2][0],
                    itemSym.getSymR()[2][1],

                    itemSym.getSymR()[2][2],
                    itemSym.getSymT()[2]);
        }
    }

    public void LATTirazation(){
        double[][] SYMr = new double[3][3];
        double[] SYMt = new double[3];
        List<SymmetryItem> subSymMass = new ArrayList<>();
        subSymMass =   (ArrayList<SymmetryItem>) deepClone(this.symMass);
        for (SymmetryItem itemSYM : subSymMass){
            switch (Math.abs(this.LATT)) {
                case 2:
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += .5;
                    SYMt[1] += .5;
                    SYMt[2] += .5;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    break;
                case 3:
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += (double) 1 / 3;
                    SYMt[1] += (double) 1 / 6;
                    SYMt[2] += (double) 1 / 6;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += (double) 1 / 6;
                    SYMt[1] += (double) 1 / 3;
                    SYMt[2] += (double) 1 / 3;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    break;
                case 4:
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += 0;
                    SYMt[1] += .5;
                    SYMt[2] += .5;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += .5;
                    SYMt[1] += 0;
                    SYMt[2] += .5;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += .5;
                    SYMt[1] += .5;
                    SYMt[2] += 0;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    break;
                case 5:
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += 0;
                    SYMt[1] += .5;
                    SYMt[2] += .5;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    break;
                case 6:
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += .5;
                    SYMt[1] += 0;
                    SYMt[2] += .5;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    break;
                case 7:
                    SYMt = (double[]) deepClone(itemSYM.getSymT());
                    SYMr = (double[][]) deepClone(itemSYM.getSymR());
                    SYMt[0] += .5;
                    SYMt[1] += .5;
                    SYMt[2] += 0;
                    this.symMass.add(new SymmetryItem(SYMr,  SYMt));
                    SYMt = null;
                    SYMr = null;
                    break;
            }
        }
    }
}

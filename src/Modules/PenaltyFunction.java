package Modules;

import CrystTools.Fragment;
import CrystTools.ReciprocalItem;
import CrystTools.Symmetry;
import CrystTools.UnitCell;
import Utilities.ObjectsUtilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by Developer on 21.09.2015.
 */
public class PenaltyFunction {
    private String penaltyDataFilename = "";

    private class parRange {
        public String OPT = "";
        public double[] Range;

        public parRange() {
            this.OPT = "";
            Range = new double[]{0.0, 0.0, 0.0, 0.0};
        }
    }

    private parRange x = new parRange();
    private parRange y = new parRange();
    private parRange z = new parRange();
    private parRange phi1 = new parRange();
    private parRange phi2 = new parRange();
    private parRange theta = new parRange();
    private parRange o = new parRange();
    private parRange u = new parRange();
    private double force = 1;
    private double speed = 1;


    public PenaltyFunction(String penaltyDataFilename) {
        this.penaltyDataFilename = penaltyDataFilename;
        List<String> input = ObjectsUtilities.getContentFromFile(this.penaltyDataFilename, "PENALTY_FUNCTION_SETTINGS");
        for (String s : input) {
            try {
                if (!s.isEmpty()) {
                    Pattern p = Pattern.compile("(\\S+)");
                    Matcher m = p.matcher(s);
                    List<String> allMatches = new ArrayList<String>();
                    while (m.find()) allMatches.add(m.group());
                    switch (allMatches.get(0)) {
                        case "FORCE":
                            this.force = Math.pow(10, Double.valueOf(allMatches.get(1)).doubleValue()) - 1;
                            break;
                        case "SPEED":
                            this.speed = Double.valueOf(allMatches.get(1)).doubleValue();
                            break;
                        case "X":
                            this.x.Range[0] = Double.valueOf(allMatches.get(1)).doubleValue();
                            this.x.Range[1] = Double.valueOf(allMatches.get(2)).doubleValue();
                            this.x.Range[2] = (this.x.Range[0] + this.x.Range[1]) / 2;
                            this.x.Range[3] = (this.x.Range[1] - this.x.Range[0]) / 2;
                            if (allMatches.size() == 4) this.x.OPT = allMatches.get(3);
                            break;
                        case "Y":
                            this.y.Range[0] = Double.valueOf(allMatches.get(1)).doubleValue();
                            this.y.Range[1] = Double.valueOf(allMatches.get(2)).doubleValue();
                            this.y.Range[2] = (this.y.Range[0] + this.y.Range[1]) / 2;
                            this.y.Range[3] = (this.y.Range[1] - this.y.Range[0]) / 2;
                            if (allMatches.size() == 4) this.y.OPT = allMatches.get(3);
                            break;
                        case "Z":
                            this.z.Range[0] = Double.valueOf(allMatches.get(1)).doubleValue();
                            this.z.Range[1] = Double.valueOf(allMatches.get(2)).doubleValue();
                            this.z.Range[2] = (this.z.Range[0] + this.z.Range[1]) / 2;
                            this.z.Range[3] = (this.z.Range[1] - this.z.Range[0]) / 2;
                            if (allMatches.size() == 4) this.z.OPT = allMatches.get(3);
                            break;
                        case "Phi1/PI":
                            this.phi1.Range[0] = Math.PI * Double.valueOf(allMatches.get(1)).doubleValue();
                            this.phi1.Range[1] = Math.PI * Double.valueOf(allMatches.get(2)).doubleValue();
                            this.phi1.Range[2] = (this.phi1.Range[0] + this.phi1.Range[1]) / 2;
                            this.phi1.Range[3] = (this.phi1.Range[1] - this.phi1.Range[0]) / 2;
                            if (allMatches.size() == 4) this.phi1.OPT = allMatches.get(3);
                            break;
                        case "Phi2/PI":
                            this.phi2.Range[0] = Math.PI * Double.valueOf(allMatches.get(1)).doubleValue();
                            this.phi2.Range[1] = Math.PI * Double.valueOf(allMatches.get(2)).doubleValue();
                            this.phi2.Range[2] = (this.phi2.Range[0] + this.phi2.Range[1]) / 2;
                            this.phi2.Range[3] = (this.phi2.Range[1] - this.phi2.Range[0]) / 2;
                            if (allMatches.size() == 4) this.phi2.OPT = allMatches.get(3);
                            break;
                        case "Theta/PI":
                            this.theta.Range[0] = Math.PI * Double.valueOf(allMatches.get(1)).doubleValue();
                            this.theta.Range[1] = Math.PI * Double.valueOf(allMatches.get(2)).doubleValue();
                            this.theta.Range[2] = (this.theta.Range[0] + this.theta.Range[1]) / 2;
                            this.theta.Range[3] = (this.theta.Range[1] - this.theta.Range[0]) / 2;
                            if (allMatches.size() == 4) this.theta.OPT = allMatches.get(3);
                            break;
                        case "O":
                            this.o.Range[0] = Double.valueOf(allMatches.get(1)).doubleValue();
                            this.o.Range[1] = Double.valueOf(allMatches.get(2)).doubleValue();
                            this.o.Range[2] = (this.o.Range[0] + this.o.Range[1]) / 2;
                            this.o.Range[3] = (this.o.Range[1] - this.o.Range[0]) / 2;
                            if (allMatches.size() == 4) this.o.OPT = allMatches.get(3);
                            break;
                        case "U":
                            this.u.Range[0] = Double.valueOf(allMatches.get(1)).doubleValue();
                            this.u.Range[1] = Double.valueOf(allMatches.get(2)).doubleValue();
                            this.u.Range[2] = (this.u.Range[0] + this.u.Range[1]) / 2;
                            this.u.Range[3] = (this.u.Range[1] - this.u.Range[0]) / 2;
                            if (allMatches.size() == 4) this.u.OPT = allMatches.get(3);
                            break;
                        default:
                            break;
                    }
                    allMatches.clear();
                }
            } catch (NumberFormatException e) {
                //throw new RuntimeException(e);
            }
        }

    }

    public double Psi1(FragmentData FRAG) {
        double sum = 0;
        for (Fragment imFragment : FRAG.getFragMass()) {
            if (!this.o.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.o.Range[3] + Math.abs(imFragment.getFragO() - this.o.Range[2]), 0), speed);
            if (!this.u.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.u.Range[3] + Math.abs(imFragment.getFragU() - this.u.Range[2]), 0), speed);
            if (!this.x.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.x.Range[3] + Math.abs(imFragment.getFragX() - this.x.Range[2]), 0), speed);
            if (!this.y.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.y.Range[3] + Math.abs(imFragment.getFragY() - this.y.Range[2]), 0), speed);
            if (!this.z.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.z.Range[3] + Math.abs(imFragment.getFragZ() - this.z.Range[2]), 0), speed);
            if (!this.phi1.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.phi1.Range[3] + Math.abs(imFragment.getFragPhi1() - this.phi1.Range[2]), 0), speed);
            if (!this.phi2.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.phi2.Range[3] + Math.abs(imFragment.getFragPhi2() - this.phi2.Range[2]), 0), speed);
            if (!this.theta.OPT.contains("S"))
                sum += Math.pow(Math.max(-this.theta.Range[3] + Math.abs(imFragment.getFragTheta() - this.theta.Range[2]), 0), speed);
        }

        return this.force * sum;
    }
}

package Modules;

import CrystTools.Fragment;
import CrystTools.Symmetry;
import CrystTools.UnitCell;
import MathTools.FastMath;
import Utilities.ObjectsUtilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static MathTools.FastMath.*;
import static Utilities.ObjectsUtilities.*;

/**
 * Created by Developer on 10.09.2015.
 * Simulated annealing module for solving strucure.
 */
public class SAmodule {
    private UnitCell CELL;
    public Symmetry SYM;
    private FragmentData FRAG;
    private DiffractionData HKL;
    private TemperatureData TEMPREGIME;
    private AnnealingSettings SASETTINGS;
    private Energy E;
    private String tempFileName;
    private int[] ParametersList;
    private String resultFileName;


    public class TemperatureData {
        private String temperatureDataFilename = "";
        private ArrayList<TemperatureItem> tempRegime = new ArrayList<>();

        public TemperatureData(String temperatureDataFilename) {
            this.temperatureDataFilename = temperatureDataFilename;
            List<String> input = ObjectsUtilities.getContentFromFile(this.temperatureDataFilename);
            for (String s : input) {
                try {
                    Pattern p = Pattern.compile("(\\S+)");
                    Matcher m = p.matcher(s);
                    List<String> allMatches = new ArrayList<String>();
                    while (m.find()) {
                        allMatches.add(m.group());
                    }
                    this.tempRegime.add(new TemperatureItem(
                            Double.valueOf(allMatches.get(0)).doubleValue(),
                            Integer.parseInt(allMatches.get(1)),
                            allMatches.get(2)
                    ));

                } catch (NumberFormatException s2nTemperatureData) {
                    //throw new RuntimeException(s2n);
                }
            }
        }

        public List<TemperatureItem> getREGIME() {
            return tempRegime;
        }

    }


    public class TemperatureItem {

        private double cycleT;
        private int cycleIterations;
        private String cycleOpt;


        public TemperatureItem(double cyclet, int cycleiterations, String cycleopt) {
            this.cycleT = cyclet;
            this.cycleIterations = cycleiterations;
            this.cycleOpt = cycleopt;

        }
    }


    public class AnnealingSettings {

        private String annealingSettingsFilename = "";
        private double P_FOR_KB = 0;
        private double RANDOMIZE_ACCURACY = 0;
        private int FIRST_RANDOMIZATION = 0;
        private double MAX_PARAMETERS_STEP = 0;
        private int PRINT_FRAGMENT = 0;
        private String SAVE_BEST_RESULT = "N";


        public AnnealingSettings(String annealingSettingsFilename) {
            this.annealingSettingsFilename = annealingSettingsFilename;
            List<String> input = ObjectsUtilities.getContentFromFile(this.annealingSettingsFilename);
            for (String s : input) {
                try {
                    if (!s.isEmpty()) {
                        Pattern p = Pattern.compile("(\\S+)");
                        Matcher m = p.matcher(s);
                        ArrayList<String> allMatches = new ArrayList<String>();
                        while (m.find()) allMatches.add(m.group());
                        switch (allMatches.get(0)) {
                            case "P_FOR_KB":
                                this.P_FOR_KB = Double.valueOf(allMatches.get(1)).doubleValue();
                                break;
                            case "RANDOMIZE_ACCURACY":
                                this.RANDOMIZE_ACCURACY = Double.valueOf(allMatches.get(1)).doubleValue();
                                break;
                            case "FIRST_RANDOMIZATION":
                                this.FIRST_RANDOMIZATION = Integer.parseInt(allMatches.get(1));
                                break;
                            case "MAX_PARAMETERS_STEP":
                                this.MAX_PARAMETERS_STEP = Double.valueOf(allMatches.get(1)).doubleValue();
                                break;
                            case "SAVE_BEST_RESULT":
                                this.SAVE_BEST_RESULT = allMatches.get(1);
                                break;
                            case "PRINT_FRAGMENT":
                                this.PRINT_FRAGMENT = Integer.parseInt(allMatches.get(1));
                                break;
                            default:
                                break;
                        }
                        allMatches.clear();
                    }
                } catch (NumberFormatException e) {
                }
            }
        }
    }


    public SAmodule(UnitCell CELL,
                    FragmentData FRAG,
                    DiffractionData HKL,
                    Energy E,
                    String annealingSettingsFilename,
                    String temperatureDataFilename,
                    String tempFileName,
                    String resultFileName) {
        this.CELL = CELL;
        this.FRAG = FRAG;
        this.HKL = HKL;
        this.TEMPREGIME = new TemperatureData(temperatureDataFilename);
        this.SASETTINGS = new AnnealingSettings(annealingSettingsFilename);
        this.tempFileName = tempFileName;
        this.resultFileName = resultFileName;
        DeleteFile(tempFileName);
        this.ParametersList = generateParametersList(FRAG);
        this.E = E;
    }

    public double calcOrder(ArrayList<Double> E, double kB, double T) {
        double mE = 0;
        double dE = 0;
        double result = 0;
        int ESIZE = E.size();
        if (ESIZE > 0) {
            for (int i = 0; i < ESIZE; i++) mE += E.get(i);
            mE /= ESIZE;
            for (int i = 0; i < ESIZE; i++) dE += Math.pow((E.get(i) - mE), 2);
            dE /= ESIZE;
            result = dE / (kB * Math.pow(T, 2));
            if (!Double.isNaN(result)) {
                return result;
            } else {
                return -1;
            }
        } else {
            return -1;
        }

    }

    public double calcBolzman(List<Double> E, double T, double springP) {
        double dE = 0;
        int ESIZE = E.size();
        if ((ESIZE - 1) > 0) {
            for (int i = 0; i < ESIZE - 1; i++) dE += E.get(i) - E.get(i + 1);
            dE /= (ESIZE - 1) / 2;
            return Math.abs(dE / (T * Math.log(springP)));
        } else {
            return 1;
        }

    }

    public double calcW(List<Double> E1, List<Double> E2) {
        double E1Mid = 0;
        double E2Mid = 0;
        if (E1.size() + E2.size() > 3) {
            for (int i = 0; i < E1.size() - 1; i++) {
                E1Mid += Math.pow(E1.get(i + 1) - E1.get(i), 2);
            }
            E1Mid /= E1.size() - 1;
            for (int i = 0; i < E2.size() - 1; i++) {
                E2Mid += Math.pow(E2.get(i + 1) - E2.get(i), 2);
            }
            E2Mid /= E2.size() - 1;
            return (E2Mid != 0) ? 0.5 * Math.sqrt(E2Mid / E1Mid) : 1;
        } else {
            return 1;
        }

    }




    public double randomizeParameters(FragmentData FRAG) {
        Random randomVal = new Random();
        int parametersChoice = randomVal.nextInt(this.ParametersList.length);
        double Delta = randomizeDouble("SYM", this.SASETTINGS.MAX_PARAMETERS_STEP);
        addToParameters(this.CELL, FRAG, this.ParametersList, parametersChoice, Delta);
        return Delta;
    }




    public String sumOptions(List<String> opt) {
        String optString = "";
        for (String itemOpt : opt) {
            optString += itemOpt;
        }
        return optString;
    }


    public void printTempInfo(FragmentData FRAG, Energy E, int loop, String regime) {
        List<String> output = new ArrayList<>();
        if (regime.equals("")) {
            output.add(String.format("%5s %4s %12s %6s %6s %6s %6s %7s %7s %7s %7s %7s %7s %7s %7s",
                    "#",
                    "TYPE",
                    "Ecore",
                    "RI",
                    "RII",
                    "RIII",
                    "RIV",
                    "X",
                    "Y",
                    "Z",
                    "Phi1",
                    "Phi2",
                    "Theta",
                    "O",
                    "U"));
        } else {
            Fragment frag = null;
            for (Fragment itemFrag : FRAG.getFragMass()) {
                if (itemFrag.getFragNum() == this.SASETTINGS.PRINT_FRAGMENT) {
                    frag = itemFrag;
                    break;
                }
            }
            output.add(String.format("%5d %4s %8e %6.3f %6.3f %6.3f %6.3f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f",
                    loop,
                    regime,
                    E.Ecore,
                    E.RI,
                    E.RII,
                    E.RIII,
                    E.RIV,
                    frag.getFragX(),
                    frag.getFragY(),
                    frag.getFragZ(),
                    frag.getFragPhi1(),
                    frag.getFragPhi2(),
                    frag.getFragTheta(),
                    frag.getFragO(),
                    frag.getFragU()));
        }
        ObjectsUtilities.putContentToFile(this.tempFileName, output, true);
    }


    public void run() {

        List<Double> statE = new ArrayList<>();
        List<Double> statEglobal = new ArrayList<>();
        List<Double> statModGsqExray = new ArrayList<>();
        List<Double> statModGsqErest = new ArrayList<>();
        List<Double> statK = new ArrayList<>();
        List<String> opt = new ArrayList<>();

        StringBuilder strOUT = new StringBuilder("");
        StringBuilder strSTATUS = new StringBuilder("");
        StringBuilder strIND = new StringBuilder("");
        double P = 0;
        double wExray = 0;
        double BOLTZMAN = 1.38064E-23;
        double kB = BOLTZMAN;
        int iterationNumber = 0;

        double deltaE = 0;
        int sizeE = 0;
        int numJumps = 0;
        int numDecreases = 0;
        int numDecreasesGlobal = 0;
        int numDecreasesLocal = 0;
        double minE = 0;
        double minEglobal = 0;
        double RI = 0;
        double RII = 0;
        double RIII = 0;

        double minEErest = 0;
        double minEExray = 0;


        strOUT.append(String.format("\n%-50s\n", "Simulated annealing method"));
        strOUT.append(String.format("%-50s\n", "+-----------------------------------------------------------------+"));
        strOUT.append(String.format("|  %-50s%12d |\n", "Number of fragments: ", FRAG.getFragMass().size()));
        strOUT.append(String.format("|  %-50s%12d |\n", "Number of used reflections: ", HKL.getHKL().size()));
        strOUT.append(String.format("|  %-50s%12d |\n", "Number of parameters: ", this.ParametersList.length));
        strOUT.append(String.format("|  %-50s%12.2e |\n", "Randomization step: ", this.SASETTINGS.MAX_PARAMETERS_STEP));
        strOUT.append(String.format("%-50s\n", "+-----------------------------------------------------------------+"));
        System.out.print(strOUT.toString());
        strOUT.setLength(0);
        System.out.println("Running...\n");

        if (this.ParametersList.length == 0) {
            E.OPT = "Xr1";
            E.calcEnergy(FRAG);
            E.printInfo();
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (S) - RI", E.RI));
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (C) - RII", E.RII));
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (W) - RIII", E.RIII));
            strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor - RIV", E.RIV));
            strOUT.append(String.format(" %-30s= %-12e\n", "Chem part of minimum E", E.Erest));
            strOUT.append(String.format(" %-30s= %-12e\n", "Energy", E.E));
            strOUT.append(String.format(" %-30s= %-12e\n", "Penalty function", E.Epenalty));
            System.out.print("\n" + strOUT + "\n");
            strOUT.setLength(0);
            strOUT.append(String.format("\n%-50s\n", "-------------------------------------------------------------------"));
            System.out.print(strOUT);
            System.out.print("Done.\n");
            FRAG.printFragsWithSym(SYM);
            System.exit(0);
        }

        if ((new File(this.resultFileName)).exists() && this.SASETTINGS.SAVE_BEST_RESULT.contains("Y")) {
            try {
                System.out.print("Fragment data loading...");
                FragmentData FRAG_LOAD = (FragmentData) recallObject(this.resultFileName);
                for (Fragment itemFrag : FRAG_LOAD.getFragMass()) {
                    itemFrag.setFragOPT(FRAG.getFragMass().get(FRAG_LOAD.getFragMass().indexOf(itemFrag)).getFragOPT());
                }
                FRAG = (FragmentData) deepClone(FRAG_LOAD);
                FRAG_LOAD = null;
            } catch (IOException | ClassNotFoundException e) {

            }
            System.out.print("\r\r");
        } else {
            randomizeParametersInitial(this.CELL, FRAG, SASETTINGS.FIRST_RANDOMIZATION);
        }
        FRAG.printFrags();
        FragmentData FRAG_I = (FragmentData) deepClone(FRAG);
        FragmentData FRAG_II = (FragmentData) deepClone(FRAG);
        FragmentData FRAG_III = (FragmentData) deepClone(FRAG);
        Iterator<TemperatureItem> iterCYCLE = this.TEMPREGIME.getREGIME().iterator();

        if (this.SASETTINGS.PRINT_FRAGMENT != 0) printTempInfo(FRAG_I, E, 0, "");

        while (iterCYCLE.hasNext()) {
            TemperatureItem itemTemperatureItem = iterCYCLE.next();

            numJumps = 0;
            numDecreases = 0;
            numDecreasesGlobal = 0;
            numDecreasesLocal = 0;
            strIND.setLength(0);
            iterationNumber++;
            statModGsqExray.clear();
            statModGsqErest.clear();
            statE.clear();
            opt.clear();

            if (itemTemperatureItem.cycleOpt.contains("Xr1")) opt.add("Xr1");
            if (itemTemperatureItem.cycleOpt.contains("Xr2")) opt.add("Xr2");
            if (itemTemperatureItem.cycleOpt.contains("Xr3")) opt.add("Xr3");
            if (itemTemperatureItem.cycleOpt.contains("Xr4")) opt.add("Xr4");
            if (itemTemperatureItem.cycleOpt.contains("Xr5")) opt.add("Xr5");
            if (itemTemperatureItem.cycleOpt.contains("Re1")) opt.add("Re1");


            FRAG_II = (FragmentData) deepClone(FRAG);
            E.OPT = sumOptions(opt);
            E.calcEnergy(FRAG_II);
            FRAG.printFrags();
            E.printInfo();

            if (!this.resultFileName.equals("") && this.SASETTINGS.SAVE_BEST_RESULT.contains("Y")) {
                try {
                    saveObject(FRAG, this.resultFileName);
                } catch (IOException e) {
                }
            }

            statE.add(E.E);
            statEglobal.add(E.E);
            statK.add(E.getK());

            minE = statE.get(statE.indexOf(Collections.min(statE)));
            minEglobal = statEglobal.get(statEglobal.indexOf(Collections.min(statEglobal)));

            if (itemTemperatureItem.cycleOpt.contains("Xr")) {
                if (itemTemperatureItem.cycleOpt.contains("K")) {
                    E.improveK(FRAG_II);
                }
            }

            if (this.SASETTINGS.PRINT_FRAGMENT != 0) printTempInfo(FRAG_II, E, 0, "*");

            for (int iteration = 0; iteration < itemTemperatureItem.cycleIterations; iteration++) {

                strSTATUS.setLength(0);
                strOUT.setLength(0);

                strOUT.append(String.format(" %-30s->%d\n", "Iteration number", iterationNumber));
                strOUT.append(String.format(" %-30s= %-12e\n", "T", itemTemperatureItem.cycleT));
                strOUT.append(String.format(" %-30s= %-12d\n", "Cycles", itemTemperatureItem.cycleIterations));


                if (itemTemperatureItem.cycleOpt.contains("S")) {
                    strSTATUS.append("Skipping loop, ");
                    break;
                }

                FRAG_I = (FragmentData) deepClone(FRAG_II);
                double deltaParameters = randomizeParameters(FRAG_I);


                if (itemTemperatureItem.cycleOpt.contains("Xr")) {
                    if (itemTemperatureItem.cycleOpt.contains("K")) {
                        strSTATUS.append("Scale factor calculation, ");
                        strOUT.append(String.format(" %-30s= %-12e\n", "Scale factor K", E.getK()));
                    }
                }


                if (itemTemperatureItem.cycleOpt.contains("Re")) {
                    if (itemTemperatureItem.cycleOpt.contains("Xr")) {
                        if (itemTemperatureItem.cycleOpt.contains("W")) {
                            wExray = calcW(statModGsqExray, statModGsqErest);
                            E.improvewExray(wExray);
                            strSTATUS.append("Weights calculation, ");
                            strOUT.append(String.format(" %-30s= %-12e\n", "Exray weight w", wExray));
                        }
                    } else {
                        wExray = 0;
                        E.improvewExray(wExray);
                        strSTATUS.append("Dynamics simulation only, ");
                    }
                }


                E.calcEnergy(FRAG_I);
                minE = statE.get(statE.indexOf(Collections.min(statE)));
                minEglobal = statEglobal.get(statEglobal.indexOf(Collections.min(statEglobal)));
                statE.add(E.E);
                statEglobal.add(E.E);



                if (deltaParameters != 0) {
                    statModGsqExray.add(Math.pow(E.Exray / deltaParameters, 2));
                    statModGsqErest.add(Math.pow(E.Erest / deltaParameters, 2));
                }


                if (itemTemperatureItem.cycleOpt.contains("Xr") || itemTemperatureItem.cycleOpt.contains("Re")) {
                    if (itemTemperatureItem.cycleOpt.contains("B")) {
                        kB = calcBolzman(statE, itemTemperatureItem.cycleT, this.SASETTINGS.P_FOR_KB);
                        strSTATUS.append("Bolzman constant optimization, ");
                        strOUT.append(String.format(" %-30s= %-12e\n", "Bolzman constant", kB));
                    }
                }



                if ((itemTemperatureItem.cycleOpt.contains("D") || itemTemperatureItem.cycleOpt.contains("J"))) {

                    if (itemTemperatureItem.cycleOpt.contains("DJ")) {
                        strSTATUS.append("SA full, ");
                    } else if (itemTemperatureItem.cycleOpt.contains("D")) {
                        strSTATUS.append("Decreases only, ");
                    } else if (itemTemperatureItem.cycleOpt.contains("J")) {
                        strSTATUS.append("Jumps only, ");
                    }

                    if (E.E < minEglobal) {
                        FRAG = (FragmentData) deepClone(FRAG_I);
                        minEExray = E.Exray;
                        minEErest = E.Erest;
                        RI = E.RI;
                        RII = E.RII;
                        RIII = E.RIII;
                        numDecreasesGlobal++;
                        strIND.setLength(0);
                        strIND.append(String.format(" %05d", numDecreasesGlobal));
                        strIND.append(String.format(" RI=%-4.2f", E.RI));
                        strIND.append(String.format(" RII=%-4.2f", E.RII));
                        strIND.append(String.format(" RIII=%-4.2f", E.RIII));
                        strIND.append(String.format(" RIV=%-4.2f", E.RIV));
                    }

                    if (E.E < minE) {
                        FRAG_III = (FragmentData) deepClone(FRAG_I);
                        numDecreasesLocal++;
                    }

                    sizeE = statE.size();
                    deltaE = statE.get(sizeE - 1) - statE.get(sizeE - 2);
                    if (deltaE < 0) {
                        if (itemTemperatureItem.cycleOpt.contains("D")) {
                            FRAG_II = (FragmentData) deepClone(FRAG_I);
                            numDecreases++;
                        }
                    } else {
                        if (itemTemperatureItem.cycleOpt.contains("J")) {
                            P = Math.exp(-deltaE / (kB * itemTemperatureItem.cycleT));
                            if (P > Math.random()) {
                                FRAG_II = (FragmentData) deepClone(FRAG_I);
                                numJumps++;
                            }
                        }
                    }

                    strSTATUS.append("Minima of E searching, ");
                    strOUT.append(String.format(" %-30s= %-12d\n", "Local jumps", numJumps));
                    strOUT.append(String.format(" %-30s= %-12d\n", "Local decreases:", numDecreasesLocal));
                    strOUT.append(String.format(" %-30s= %-12d\n", "Global decreases:", numDecreasesGlobal));
                    strOUT.append(String.format(" %-30s= %-12e\n", "Mininmum E", minE));
                    if (RI != 0) strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (with scale)", RI));
                    if (RII != 0) strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (without scale)", RII));
                    if (RIII != 0) strOUT.append(String.format(" %-30s= %-12.3f\n", "R-factor (with weights)", RIII));
                    if (minEExray != 0)
                        strOUT.append(String.format(" %-30s= %-12e\n", "Xray part of minimum E", minEExray));
                    if (minEErest != 0)
                        strOUT.append(String.format(" %-30s= %-12e\n", "Chem part of minimum E", minEErest));

                }
                System.out.print("\r" +
                        strSTATUS.substring(0, strSTATUS.length() - 2) +
                        ": " +
                        String.format("%2.0f", Double.valueOf(100 * (iteration + 1) / itemTemperatureItem.cycleIterations)) +
                        "% " + strIND);
            }

            System.out.print("\r");
            System.out.print("\n" + strOUT + "\n");
        }

        strOUT.append(String.format("%-50s\n", "-------------------------------------------------------------------"));
        System.out.print(strOUT);
        System.out.print("Done\n");
    }
}

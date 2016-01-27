package Modules;

import CrystTools.*;
import MathTools.FastMath;
import Utilities.ObjectsUtilities;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static Utilities.ObjectsUtilities.generateAtomNum;
import static Utilities.ObjectsUtilities.readFileToList;

/**
 * Created by Developer on 10.09.2015.
 */
public class FragmentData implements Serializable {
    private String fragDataTypesFilename = "";
    private String fragDataListFilename = "";
    private String fragDataPrintListFilename = "";
    private String aDataScatteringFilename = "";
    private String aDataAtomParamsFilename = "";
    private UnitCell CELL;
    private List<Fragment> fragMass = new ArrayList<>();

    public FragmentData(UnitCell CELL,
                        String fragDataTypesFilename,
                        String fragDataListFilename,
                        String fragDataPrintListFilename,
                        String aDataScatteringFilename,
                        String aDataAtomParamsFilename) {
        this.CELL = CELL;
        this.fragDataTypesFilename = fragDataTypesFilename;
        this.fragDataListFilename = fragDataListFilename;
        this.fragDataPrintListFilename = fragDataPrintListFilename;
        this.aDataScatteringFilename = aDataScatteringFilename;
        this.aDataAtomParamsFilename = aDataAtomParamsFilename;

        List<String> input = ObjectsUtilities.getContentFromFile(fragDataListFilename);
        for (String s : input) {
            try {
                Pattern p = Pattern.compile("(\\S+)");
                Matcher m = p.matcher(s);
                List<String> allMatches = new ArrayList<String>();
                while (m.find()) {
                    allMatches.add(m.group());
                }
                String fName = allMatches.get(0);
                List<Atom> fAtoms = findAtomsCoords(fName);
                int fNum = Integer.valueOf(allMatches.get(1));
                double fO = Double.valueOf(allMatches.get(2)).doubleValue();
                double fU = Double.valueOf(allMatches.get(3)).doubleValue();
                String fSt = allMatches.get(4);
                String fOpt = "";
                if (allMatches.size() == 6) fOpt = allMatches.get(5);
                if (!fOpt.contains("S")) this.fragMass.add(new Fragment(fName, fAtoms, fNum, fO, fU, fSt, fOpt));
            } catch (NumberFormatException e) {
                throw new RuntimeException(e);
            }
        }
    }

    public List<Atom> findAtomsCoords(String fragmentName) {
        List<Atom> fragAtoms = new ArrayList<>();
        List<String> input = ObjectsUtilities.getContentFromFile(this.fragDataTypesFilename);
        boolean flag = false;
        for (String s : input) {
            if (s.equals(fragmentName)) {
                flag = true;
                continue;
            }
            if (s.equals("END")) flag = false;
            if (flag) {
                try {
                    Pattern p = Pattern.compile("(\\S+)");
                    Matcher m = p.matcher(s);
                    List<String> allMatches = new ArrayList<>();
                    while (m.find()) {
                        allMatches.add(m.group());
                    }
                    String aName = allMatches.get(0);
                    double[] VECT = this.CELL.c2f(
                            Double.valueOf(allMatches.get(1)).doubleValue(),
                            Double.valueOf(allMatches.get(2)).doubleValue(),
                            Double.valueOf(allMatches.get(3)).doubleValue()
                    );
                    double[] aScat = findAtomsScattering(aName);
                    double[] aParams = findAtomsParams(aName);
                    fragAtoms.add(new Atom(aName,
                            VECT[0],
                            VECT[1],
                            VECT[2],
                            aScat,
                            aParams[0],
                            aParams[1]));
                } catch (NumberFormatException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        return fragAtoms;
    }

    public double[] findAtomsScattering(String aName) {
        double[] aScat = new double[9];
        List<String> input = ObjectsUtilities.getContentFromFile(this.aDataScatteringFilename);
        for (String s : input) {
            try {
                Pattern p = Pattern.compile("(\\S+)");
                Matcher m = p.matcher(s);
                List<String> allMatches = new ArrayList<>();
                while (m.find()) {
                    allMatches.add(m.group());
                }
                if (allMatches.get(0).equals(aName)) {
                    for (int i = 0; i < 9; i++) aScat[i] = Double.valueOf(allMatches.get(i + 1)).doubleValue();
                }
            } catch (NumberFormatException e) {
                throw new RuntimeException(e);
            }
        }
        return aScat;
    }

    public double[] findAtomsParams(String aName) {
        double[] aParams = new double[2];
        List<String> input = ObjectsUtilities.getContentFromFile(this.aDataAtomParamsFilename);
        for (String s : input) {
            try {
                Pattern p = Pattern.compile("(\\S+)");
                Matcher m = p.matcher(s);
                ArrayList<String> allMatches = new ArrayList<>();
                while (m.find()) {
                    allMatches.add(m.group());
                }
                if (allMatches.get(0).equals(aName)) {
                    aParams[0] = Double.valueOf(allMatches.get(1)).doubleValue();
                    aParams[1] = Double.valueOf(allMatches.get(2)).doubleValue();
                }
            } catch (NumberFormatException e) {
                throw new RuntimeException(e);
            }
        }
        return aParams;
    }


    public List<Fragment> getFragMass() {
        return fragMass;
    }


    public void printFrags() {
        printFrags(this.fragDataPrintListFilename);
    }


    public void printFrags(String fragDataPrintListFilename) {
        int numC = 0;
        int numM = 0;
        List<String> output = new ArrayList<>();
        List<String> aTypeCounter = new ArrayList<>();
        for (Fragment itemFrag : this.fragMass) {
            ++numM;
            output.add(String.format("MOLE %d", numM));
            output.add(String.format("AFIX 1"));
            for (Atom itemAtom : itemFrag.getFragAtoms()) {
                ++numC;
                if (!aTypeCounter.contains(itemAtom.getAtomName())) aTypeCounter.add(itemAtom.getAtomName());
                output.add(String.format("%-4s  %d  % .5f  % .5f  % .5f 1%.4f %.3f",
                        generateAtomNum(itemAtom.getAtomName(), numC),
                        aTypeCounter.indexOf(itemAtom.getAtomName()) + 1,
                        itemAtom.getAtomX(),
                        itemAtom.getAtomY(),
                        itemAtom.getAtomZ(),
                        Math.abs(itemFrag.getFragO()),
                        itemFrag.getFragU()));
            }
            output.add(String.format("AFIX 0"));
        }

        ObjectsUtilities.putContentToFile(fragDataPrintListFilename, output);

    }

    public void printFragsWithSym(Symmetry sym) {
        printFragsWithSym(this.fragDataPrintListFilename, sym);
    }


    public void printFragsWithSym(String fragDataPrintListFilename, Symmetry SYM) {
        int numC = 0;
        int numM = 0;
        List<String> output = new ArrayList<>();
        List<String> aTypeCounter = new ArrayList<>();
        for (Fragment itemFrag : this.fragMass) {
            for (SymmetryItem itemSYM : SYM.getSymMass()) {
                ++numM;
                output.add(String.format("MOLE %d", numM));
                output.add(String.format("AFIX 1"));
                for (Atom itemAtom : itemFrag.getFragAtoms()) {
                    ++numC;
                    double[] VECTCORD = {itemAtom.getAtomX(), itemAtom.getAtomY(), itemAtom.getAtomZ()};
                    double[] VECT = FastMath.VpV(itemSYM.getSymT(), FastMath.MmV(itemSYM.getSymR(), VECTCORD));
                    if (!aTypeCounter.contains(itemAtom.getAtomName()))
                        aTypeCounter.add(itemAtom.getAtomName());
                    output.add(String.format("%-4s  %d  % .5f  % .5f  % .5f 1%.4f %.3f",
                            generateAtomNum(itemAtom.getAtomName(), numC),
                            aTypeCounter.indexOf(itemAtom.getAtomName()) + 1,
                            VECT[0],
                            VECT[1],
                            VECT[2],
                            Math.abs(itemFrag.getFragO()),
                            itemFrag.getFragU()));
                }
                output.add(String.format("AFIX 0"));
            }
        }
        ObjectsUtilities.putContentToFile(fragDataPrintListFilename, output);

    }


    public void printFragParameters(String fragParametersPrintListFilename, int fragNum, String additionalString) {
        List<String> output = new ArrayList<>();
        int indexFrag = 0;
        for (Fragment itemFrag : this.fragMass) {
            if (itemFrag.getFragNum() == fragNum) {
                indexFrag = this.fragMass.indexOf(itemFrag);
                break;
            }
        }
        output.add(String.format("%s % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f",
                additionalString,
                this.getFragMass().get(indexFrag).getFragX(),
                this.getFragMass().get(indexFrag).getFragY(),
                this.getFragMass().get(indexFrag).getFragZ(),
                this.getFragMass().get(indexFrag).getFragPhi1(),
                this.getFragMass().get(indexFrag).getFragPhi2(),
                this.getFragMass().get(indexFrag).getFragTheta(),
                this.getFragMass().get(indexFrag).getFragO(),
                this.getFragMass().get(indexFrag).getFragU()));
        ObjectsUtilities.putContentToFile(fragParametersPrintListFilename, output, true);
    }

    public void fragsParametersAdjustment() {
        for (Fragment itemFrag : this.fragMass) {
            itemFrag.fragParametersToOrder();
        }
    }

}

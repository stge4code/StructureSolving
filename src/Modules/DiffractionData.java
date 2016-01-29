package Modules;

import CrystTools.Atom;
import CrystTools.Fragment;
import CrystTools.ReciprocalItem;
import CrystTools.UnitCell;
import MathTools.FastMath;
import Utilities.ObjectsUtilities;

import java.io.*;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static Utilities.ObjectsUtilities.generateAtomNum;

/**
 * Created by Developer on 10.09.2015.
 */
public class DiffractionData {
    private String diffDataFilename = "";
    private String diffractionDataFilenameExport = "";
    private List<ReciprocalItem> HKL = new ArrayList<>();
    private DiffractionDataSettings DIFDATASETTINGS;


    public class DiffractionDataSettings {
        private String diffractionDataSettingsFilename;
        private String MERGE;
        private String ORDER;
        private double SCATTERING_MERGE_ACCURACY;
        private int[] N = new int[2];
        private SortRules SORT_H = new SortRules();
        private SortRules SORT_K = new SortRules();
        private SortRules SORT_L = new SortRules();
        private SortRules SORT_RES = new SortRules();
        private SortRules SORT_ISI = new SortRules();
        private SortRules SORT_I = new SortRules();

        private class SortRules {
            public ArrayList<String> l = new ArrayList<>();
            public ArrayList<String> e = new ArrayList<>();
            public ArrayList<String> g = new ArrayList<>();
        }

        private SortRules genSortRules(String S) {
            SortRules SORT = new SortRules();
            Matcher mS = Pattern.compile("[^,]+").matcher(S);
            ArrayList<String> allMatchesS = new ArrayList<String>();
            while (mS.find()) allMatchesS.add(mS.group());
            for (String itemS : allMatchesS) {
                if (itemS.substring(0, 1).equals(">")) SORT.g.add(itemS.substring(1));
                if (itemS.substring(0, 1).equals("=")) SORT.e.add(itemS.substring(1));
                if (itemS.substring(0, 1).equals("<")) SORT.l.add(itemS.substring(1));
            }
            return SORT;
        }

        public DiffractionDataSettings(String diffractionDataFilename) {
            this.diffractionDataSettingsFilename = diffractionDataFilename;
            List<String> input = ObjectsUtilities.getContentFromFile(diffractionDataFilename);
            for (String s : input) {
                try {
                    if (!s.isEmpty()) {
                        Pattern p = Pattern.compile("(\\S+)");
                        Matcher m = p.matcher(s);
                        List<String> allMatches = new ArrayList<String>();
                        while (m.find()) allMatches.add(m.group());
                        switch (allMatches.get(0)) {
                            case "RES":
                                this.SORT_RES = genSortRules(allMatches.get(1));
                                break;
                            case "SCATTERING_MERGE_ACCURACY":
                                this.SCATTERING_MERGE_ACCURACY = Double.valueOf(allMatches.get(1)).doubleValue();
                                break;
                            case "I/S(I)":
                                this.SORT_ISI = genSortRules(allMatches.get(1));
                                break;
                            case "I":
                                this.SORT_I = genSortRules(allMatches.get(1));
                                break;
                            case "H":
                                this.SORT_H = genSortRules(allMatches.get(1));
                                break;
                            case "K":
                                this.SORT_K = genSortRules(allMatches.get(1));
                                break;
                            case "L":
                                this.SORT_L = genSortRules(allMatches.get(1));
                                break;
                            case "N":
                                this.N[0] = (int) Double.valueOf(allMatches.get(1)).doubleValue();
                                this.N[1] = (int) Double.valueOf(allMatches.get(2)).doubleValue();
                                break;
                            case "MERGE":
                                this.MERGE = allMatches.get(1);
                                break;
                            case "ORDER":
                                this.ORDER = allMatches.get(1);
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

    public DiffractionData(UnitCell CELL,
                           String diffdatafilename,
                           String diffractionSettingsFilename,
                           String diffractionDataFilenameExport) {

        this.diffDataFilename = diffdatafilename;
        this.DIFDATASETTINGS = new DiffractionDataSettings(diffractionSettingsFilename);
        this.diffractionDataFilenameExport = diffractionDataFilenameExport;

        List<String> input = ObjectsUtilities.getContentFromFile(diffDataFilename);
        for (String s : input) {
            try {
                Matcher m = Pattern.compile("^(.{4})(.{4})(.{4})(.{8})(.{8})(.*)$").matcher(s);
                if (m.find()) {
                    if (m.group(4) != null) {
                        int h = (int) Double.valueOf(m.group(1)).doubleValue();
                        int k = (int) Double.valueOf(m.group(2)).doubleValue();
                        int l = (int) Double.valueOf(m.group(3)).doubleValue();
                        double Fsq = Double.valueOf(m.group(4)).doubleValue();
                        double sigmaFsq = Double.valueOf(m.group(5)).doubleValue();
                        double batchNumber = (!m.group(6).isEmpty()) ? Double.valueOf(m.group(6)).doubleValue() : 1.0;
                        double scatvect = CELL.calcScatVect(h, k, l);
                        this.HKL.add(new ReciprocalItem(h, k, l, Fsq, sigmaFsq, scatvect, batchNumber));
                    }
                }
            } catch (NumberFormatException e) {
            }
        }

        System.out.print("\rLoading reflections...");
        this.HKL = sortReflections(CELL, this.HKL);

        if (!DIFDATASETTINGS.MERGE.equals("N")) {
            System.out.print("\rMerging reflections...");
            for (char mch : DIFDATASETTINGS.MERGE.toCharArray()) {
                switch (mch) {
                    case 'F' :
                        this.HKL = mergeFriedelReflections(this.HKL);
                        break;
                    case 'U' :
                        this.HKL = mergeUniqueReflections(this.HKL);
                        break;
                    case 'S' :
                        this.HKL = mergeSameScatteringReflections(this.HKL);
                        break;
                   default:
                       break;
                }
            }
        }

        if (!DIFDATASETTINGS.ORDER.equals("N")) {
            System.out.print("\rOrdering reflections...");
            for (char mch : DIFDATASETTINGS.ORDER.toCharArray()) {
                switch (mch) {
                    case 'D' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.scatvect > b.scatvect)? 1 : ((a.scatvect == b.scatvect) ? 0 : -1));
                                });
                        break;
                    case 'd' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.scatvect > b.scatvect)? -1 : ((a.scatvect == b.scatvect) ? 0 : 1));
                                });
                        break;
                    case 'H' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.h > b.h)? -1 : ((a.h == b.h) ? 0 : 1));
                                });
                        break;
                    case 'h' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.h > b.h)? 1 : ((a.h == b.h) ? 0 : -1));
                                });
                        break;
                    case 'L' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.l > b.l)? -1 : ((a.l == b.l) ? 0 : 1));
                                });
                        break;
                    case 'l' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.l > b.l)? 1 : ((a.l == b.l) ? 0 : -1));
                                });
                        break;
                    case 'K' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.k > b.k)? -1 : ((a.k == b.k) ? 0 : 1));
                                });
                        break;
                    case 'k' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.k > b.k)? 1 : ((a.k == b.k) ? 0 : -1));
                                });
                        break;
                    case 'I' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.Fsq > b.Fsq)? -1 : ((a.Fsq == b.Fsq) ? 0 : 1));
                                });
                        break;
                    case 'i' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.Fsq > b.Fsq)? 1 : ((a.Fsq == b.Fsq) ? 0 : -1));
                                });
                        break;
                    case 'S' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.sigmaFsq > b.sigmaFsq)? -1 : ((a.sigmaFsq == b.sigmaFsq) ? 0 : 1));
                                });
                        break;
                    case 's' :
                        Collections.sort(this.HKL,
                                (ReciprocalItem a, ReciprocalItem b) -> {
                                    return ((a.sigmaFsq > b.sigmaFsq)? 1 : ((a.sigmaFsq == b.sigmaFsq) ? 0 : -1));
                                });
                        break;
                    default:
                        break;
                }
            }
        }


        System.out.print("\r\r");

        if (this.HKL.isEmpty()) System.exit(0);
    }


    private List<ReciprocalItem> mergeUniqueReflections(List<ReciprocalItem> unMerged) {
        List<ReciprocalItem> merged = new ArrayList<>();
        List<Integer> counter = new ArrayList<>();
        for (Iterator<ReciprocalItem> iterHKL = unMerged.iterator(); iterHKL.hasNext(); ) {
            ReciprocalItem itemHKL = iterHKL.next();
            int i = 0;
            for (ReciprocalItem itemHKLm : merged) {
                if ((itemHKLm.h == itemHKL.h) && (itemHKLm.k == itemHKL.k) && (itemHKLm.l == itemHKL.l)) {
                    int i_ = counter.get(i).intValue();
                    counter.set(i, Integer.valueOf(i_ + 1));
                    itemHKLm.Fsq *= i_;
                    itemHKLm.Fsq += itemHKL.Fsq;
                    itemHKLm.Fsq /= (i_ + 1.0);
                    itemHKLm.sigmaFsq *= i_;
                    itemHKLm.sigmaFsq += itemHKL.sigmaFsq;
                    itemHKLm.sigmaFsq /= (i_ + 1.0);
                    break;
                }
                i++;
            }
            if (i == merged.size()) {
                merged.add(new ReciprocalItem(itemHKL));
                counter.add(Integer.valueOf(1));
            }
        }
        return merged;
    }

    private List<ReciprocalItem> mergeFriedelReflections(List<ReciprocalItem> unMerged) {
        List<ReciprocalItem> merged = new ArrayList<>();
        List<Integer> counter = new ArrayList<>();
        for (Iterator<ReciprocalItem> iterHKL = unMerged.iterator(); iterHKL.hasNext(); ) {
            ReciprocalItem itemHKL = iterHKL.next();
            int i = 0;
            for (ReciprocalItem itemHKLm : merged) {
                if ((itemHKLm.h == -itemHKL.h) && (itemHKLm.k == -itemHKL.k) && (itemHKLm.l == -itemHKL.l)) {
                    int i_ = counter.get(i).intValue();
                    counter.set(i, Integer.valueOf(i_ + 1));
                    itemHKLm.Fsq *= i_;
                    itemHKLm.Fsq += itemHKL.Fsq;
                    itemHKLm.Fsq /= (i_ + 1.0);
                    itemHKLm.sigmaFsq *= i_;
                    itemHKLm.sigmaFsq += itemHKL.sigmaFsq;
                    itemHKLm.sigmaFsq /= (i_ + 1.0);
                    break;
                }
                i++;
            }
            if (i == merged.size()) {
                merged.add(new ReciprocalItem(itemHKL));
                counter.add(Integer.valueOf(1));
            }
        }
        return merged;
    }


    private List<ReciprocalItem> mergeSameScatteringReflections(List<ReciprocalItem> unMerged) {
        List<ReciprocalItem> merged = new ArrayList<>();
        List<Integer> counter = new ArrayList<>();
        for (Iterator<ReciprocalItem> iterHKL = unMerged.iterator(); iterHKL.hasNext(); ) {
            ReciprocalItem itemHKL = iterHKL.next();
            int i = 0;
            for (ReciprocalItem itemHKLm : merged) {
                if (Math.abs(itemHKLm.scatvect - itemHKL.scatvect) < this.DIFDATASETTINGS.SCATTERING_MERGE_ACCURACY) {
                    if ((itemHKL.h <= itemHKLm.h) && (itemHKL.k <= itemHKLm.k) && (itemHKL.l <= itemHKLm.l)) {
                        itemHKLm.h = itemHKL.h;
                        itemHKLm.k = itemHKL.k;
                        itemHKLm.l = itemHKL.l;
                    }
                    int i_ = counter.get(i).intValue();
                    counter.set(i, Integer.valueOf(i_ + 1));
                    itemHKLm.Fsq += itemHKL.Fsq;
                    itemHKLm.sigmaFsq *= i_;
                    itemHKLm.sigmaFsq += itemHKL.sigmaFsq;
                    itemHKLm.sigmaFsq /= (i_ + 1.0);
                    break;
                }
                i++;
            }
            if (i == merged.size()) {
                merged.add(new ReciprocalItem(itemHKL));
                counter.add(Integer.valueOf(1));
            }
        }
        return merged;
    }



    private ArrayList<ReciprocalItem> sortReflections(UnitCell CELL, List<ReciprocalItem> unSorted) {
        ArrayList<ReciprocalItem> sorted = new ArrayList<>();
        for (ReciprocalItem itemHKL : unSorted) {
            if (checkRules(CELL, itemHKL)) sorted.add(new ReciprocalItem(itemHKL));
        }
        return sorted;
    }

    private boolean checkRules(UnitCell CELL, ReciprocalItem HKL) {
        boolean condition = true;
        double checkvalue = 0;
        double resolution = 0;
        try {

            for (String S : DIFDATASETTINGS.SORT_I.l) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.Fsq < checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_I.e) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.Fsq == checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_I.g) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.Fsq  > checkvalue);
                if (!condition) return false;
            }

            for (String S : DIFDATASETTINGS.SORT_ISI.l) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.Fsq / HKL.sigmaFsq < checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_ISI.e) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.Fsq / HKL.sigmaFsq == checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_ISI.g) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.Fsq / HKL.sigmaFsq > checkvalue);
                if (!condition) return false;
            }
            resolution = 1.0 / 2.0 / HKL.scatvect;

            for (String S : DIFDATASETTINGS.SORT_RES.l) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (resolution < checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_RES.e) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (resolution == checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_RES.g) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (resolution > checkvalue);
                if (!condition) return false;
            }

            for (String S : DIFDATASETTINGS.SORT_H.l) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.h < checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_H.e) {
                if (S.contains("N")) {
                    boolean condition_temp = false;
                    for (int N = this.DIFDATASETTINGS.N[0]; N < this.DIFDATASETTINGS.N[1] - 1; N++) {
                        condition_temp = condition_temp | (FastMath.eval((S)
                                .replaceAll("N", "(" + Integer.toString(N)) + ")") - HKL.h == 0);
                        if (condition_temp) break;
                    }
                    condition = condition & condition_temp;
                } else if (S.contains("L") | S.contains("K")) {
                    condition = condition & (FastMath.eval((S)
                            .replaceAll("L", "(" + Integer.toString(HKL.l) + ")")
                            .replaceAll("K", "(" + Integer.toString(HKL.k) + ")")) - HKL.h == 0);
                } else {
                    checkvalue = Double.valueOf(S).doubleValue();
                    condition = condition & (HKL.h == checkvalue);
                    if (!condition) return false;
                }
            }
            for (String S : DIFDATASETTINGS.SORT_H.g) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.h > checkvalue);
                if (!condition) return false;

            }

            for (String S : DIFDATASETTINGS.SORT_K.l) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.k < checkvalue);
                if (!condition) return false;

            }
            for (String S : DIFDATASETTINGS.SORT_K.e) {
                if (S.contains("N")) {
                    boolean condition_temp = false;
                    for (int N = this.DIFDATASETTINGS.N[0]; N < this.DIFDATASETTINGS.N[1] - 1; N++) {
                        condition_temp = condition_temp | (FastMath.eval((S)
                                .replaceAll("N", "(" + Integer.toString(N) + ")")) - HKL.k == 0);
                        if (condition_temp) break;
                    }
                    condition = condition & condition_temp;
                } else if (S.contains("L") | S.contains("H")) {
                    condition = condition & (FastMath.eval((S)
                            .replaceAll("H", "(" + Integer.toString(HKL.h) + ")")
                            .replaceAll("L", "(" + Integer.toString(HKL.l) + ")")) - HKL.k == 0);
                } else {
                    checkvalue = Double.valueOf(S).doubleValue();
                    condition = condition & (HKL.k == checkvalue);
                    if (!condition) return false;
                }
            }
            for (String S : DIFDATASETTINGS.SORT_K.g) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.k > checkvalue);
                if (!condition) return false;
            }


            for (String S : DIFDATASETTINGS.SORT_L.l) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.l < checkvalue);
                if (!condition) return false;
            }
            for (String S : DIFDATASETTINGS.SORT_L.e) {
                if (S.contains("N")) {
                    boolean condition_temp = false;
                    for (int N = this.DIFDATASETTINGS.N[0]; N < this.DIFDATASETTINGS.N[1] - 1; N++) {
                        condition_temp = condition_temp | (FastMath.eval((S)
                                .replaceAll("N", "(" + Integer.toString(N) + ")")) - HKL.l == 0);
                        if (condition_temp) break;
                    }
                    condition = condition & condition_temp;
                } else if (S.contains("H") | S.contains("K")) {
                    condition = condition & (FastMath.eval((S)
                            .replaceAll("H", "(" + Integer.toString(HKL.h) + ")")
                            .replaceAll("K", "(" + Integer.toString(HKL.k) + ")")) - HKL.l == 0);
                } else {
                    checkvalue = Double.valueOf(S).doubleValue();
                    condition = condition & (HKL.l == checkvalue);
                    if (!condition) return false;
                }

            }
            for (String S : DIFDATASETTINGS.SORT_L.g) {
                checkvalue = Double.valueOf(S).doubleValue();
                condition = condition & (HKL.l > checkvalue);
                if (!condition) return false;
            }
        } catch (NumberFormatException e) {
            System.out.println("Incorrect input files...");
            System.exit(0);
        }
        return condition;
    }

    public List<ReciprocalItem> getHKL() {
        return HKL;
    }

    public void printHKL() {

        List<String> output = new ArrayList<>();
        for (ReciprocalItem itemHKL : this.HKL) {
            output.add(String.format("% 4d% 4d% 4d %8g %8g% 4d",
                    itemHKL.h,
                    itemHKL.k,
                    itemHKL.l,
                    FastMath.round(itemHKL.Fsq, 6),
                    FastMath.round(itemHKL.sigmaFsq, 6),
                    (int) itemHKL.batchNumber));
        }
        ObjectsUtilities.putContentToFile(this.diffractionDataFilenameExport, output);
    }
}

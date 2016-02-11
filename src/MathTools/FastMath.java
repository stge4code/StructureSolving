package MathTools;

import CrystTools.Fragment;
import CrystTools.ReciprocalItem;
import CrystTools.UnitCell;
import Modules.FragmentData;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * Created by Developer on 10.09.2015.
 */
public final class FastMath {

    private FastMath() {
    }

    public static double d2rSin(double degree) {
        return Math.sin(Math.toRadians(degree));
    }

    public static double d2rCos(double degree) {
        return Math.cos(Math.toRadians(degree));
    }

    public static double d2rTg(double degree) {
        return Math.tan(Math.toRadians(degree));
    }

    public static double[] MmV(double[][] M,
                               double[] V) {
        double[] VECT = new double[3];
        VECT[0] = M[0][0] * V[0] + M[0][1] * V[1] + M[0][2] * V[2];
        VECT[1] = M[1][0] * V[0] + M[1][1] * V[1] + M[1][2] * V[2];
        VECT[2] = M[2][0] * V[0] + M[2][1] * V[1] + M[2][2] * V[2];
        return VECT;
    }


    public static double[][] MmM(double[][] M1, double[][] M2) {
        double[][] M3 = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    M3[i][j] += M1[i][k] * M2[k][j];
                }
            }
        }
        return M3;
    }

    public static double[] VpV(double[] V1, double[] V2) {
        double[] V3 = new double[3];
        for (int i = 0; i < 3; i++) {
            V3[i] = V1[i] + V2[i];
        }
        return V3;
    }

    public static boolean allDifferent(double[] parameters) {
        boolean result = true;
        for (int i = 0; i < parameters.length; i++)
            for (int k = 0; k < parameters.length; k++) {
                result = result & ((parameters[i] != parameters[k]) | (i == k));
            }
        return result;
    }

    public static boolean allEqual(double[] parameters) {
        boolean result = true;
        for (int i = 1; i < parameters.length; i++) result = result & (parameters[i - 1] == parameters[i]);
        return result;
    }

    public static boolean AnyEquals(ReciprocalItem[] parameters) {
        boolean result = false;
        for (int i = 0; i < parameters.length; i++)
            for (int k = 0; k < parameters.length; k++) {
                if (i != k) result = result | ReciprocalItem.compare(parameters[i], parameters[k]);

            }
        return result;
    }

    public static boolean AnyEquals(ReciprocalItem parameter, ReciprocalItem[] parameters) {
        boolean result = false;
        for (int i = 0; i < parameters.length; i++) {
            result = result | ReciprocalItem.compare(parameters[i], parameter);
        }
        return result;
    }

    public static double[] KmV(double K, double[] V1) {
        double[] V2 = new double[V1.length];
        for (int i = 0; i < V1.length; i++) V2[i] = K * V1[i];
        return V2;
    }

    public static double modV(double[] V) {
        double mV = 0;
        for (int i = 0; i < V.length; i++) {
            mV += Math.pow(V[i], 2);
        }
        return Math.sqrt(mV);
    }

    public static double VmV(double[] V1, double[] V2) {
        double V3 = 0;
        for (int i = 0; i < 3; i++) {
            V3 += V1[i] * V2[i];
        }
        return V3;
    }


    public static double round(double value, int places) {
        if (Double.isNaN(value)) return Double.NaN;
        if (places < 0) throw new IllegalArgumentException();
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }


    public static int minParInV(double[] V) {
        double min = V[0];
        for (int i = 0; i < V.length; i++) min = Math.min(min, V[i]);
        for (int i = 0; i < V.length; i++) if (min == V[i]) return i;
        return -1;
    }
//    public static ArrayList<BigInteger> findFibonacciNumbers(int N){
//        ArrayList<BigInteger> F = new ArrayList<>();
//        BigInteger middleFibb = BigInteger.valueOf(N);
//        F.add(BigInteger.ONE);
//        F.add(BigInteger.ONE);
//        int i = 1;
//        do {
//            F.add(F.get(i-1).add(F.get(i)));
//            i++;
//        } while ( (F.get(i).compareTo(middleFibb) == -1 ) || ( F.get(i-1).compareTo(middleFibb) == 1) );
//        return  F;
//    }


//    public static ArrayList<BigInteger> findFibonacciNumbers(int N){
//        ArrayList<BigInteger> F = new ArrayList<>();
//        BigInteger middleFibb = BigInteger.valueOf(N);
//        F.add(BigInteger.ONE);
//        F.add(BigInteger.ONE);
//        for ( int i = 1; i < N; i++ ) F.add(F.get(i-1).add(F.get(i)));
//        return  F;
//    }
//
//    public static double findFibonacciNumbersRation(BigInteger A, BigInteger B){
//        BigDecimal Ad = new BigDecimal(A);
//        BigDecimal Bd = new BigDecimal(B);
//        Ad.divide(Bd);
//        return  Ad.doubleValue();
//    }


    public static ArrayList<Double> findFibonacciNumbersRatios(int N) {
        ArrayList<BigDecimal> F = new ArrayList<>();
        ArrayList<Double> F_Ratios = new ArrayList<>();
        int MAX_F = Math.max(3, N);
        F.add(BigDecimal.ONE);
        F.add(BigDecimal.ONE);
        for (int i = 1; i < MAX_F; i++) F.add(F.get(i - 1).add(F.get(i)));
        for (int i = 0; i < MAX_F; i++)
            F_Ratios.add((F.get(i).divide(F.get(MAX_F - 1), new MathContext(25, RoundingMode.FLOOR))).doubleValue());
        return F_Ratios;
    }


    public static Long findFibonacciNumber(int n) {
        if (n < 2)
            return 1L;
        return findFibonacciNumber(n - 1) + findFibonacciNumber(n - 2);
    }


    public static double returnInRange(double PAR, double PERIOD) {
        if (PAR > PERIOD) return (PAR % PERIOD);
        if (PAR < 0) return (PERIOD - Math.abs(PAR) % PERIOD);
        return PAR;
    }

    public static double accordingToRange(int mode, double PAR, double dPAR, double limit_left, double limit_right) {
        if (mode == 0) {
            if (PAR + dPAR > limit_right) return limit_right;
            if (PAR + dPAR < limit_left) return limit_left;
            return PAR + dPAR;
        } else if (mode == 1) {
            if (PAR + dPAR > limit_right) return limit_right;
            return PAR + dPAR;
        } else if (mode == -1) {
            if (PAR + dPAR < limit_left) return limit_left;
            return PAR + dPAR;
        } else {
            return PAR;
        }
    }


    public static double parameterRangeReflection(double PAR, double[] RANGE) {
        if (PAR > RANGE[1]) return 2 * RANGE[1] - PAR;
        if (PAR < RANGE[0]) return 2 * RANGE[0] - PAR;
        return PAR;
    }


    public static double[] c2s(double[] VECT) {
        return c2s(VECT[0], VECT[1], VECT[2]);
    }

    public static double[] c2s(double X, double Y, double Z) {
        double[] VECT = new double[3];
        VECT[0] = Math.sqrt(Math.pow(X, 2) + Math.pow(Y, 2) + Math.pow(Z, 2));
        VECT[1] = Math.acos(Z / VECT[0]);

        if (X > 0) {
            VECT[2] = Math.atan(Y / X);
        } else if (X < 0) {
            if (Y > 0) {
                VECT[2] = Math.PI + Math.atan(Y / X);
            } else if (Y < 0) {
                VECT[2] = -(Math.PI - Math.atan(Y / X));
            } else if (Y == 0) {
                VECT[2] = Math.PI;
            }
        } else if (X == 0) {
            if (Y > 0) {
                VECT[2] = Math.PI / 2;
            } else if (Y < 0) {
                VECT[2] = -Math.PI / 2;
            } else if (Y == 0) {
                VECT[2] = 0;
            }
        }
        return VECT;
    }

    public static double[] s2c(double[] VECT) {
        return s2c(VECT[0], VECT[1], VECT[2]);
    }

    public static double[] s2c(double R, double T, double P) {
        double[] VECT = new double[3];
        if (T > Math.PI) T = 2 * Math.PI - T;
        if (T < 0) T = -T;
        if (P > Math.PI) P = P - 2 * Math.PI;
        if (P < -Math.PI) P = 2 * Math.PI + P;
        VECT[0] = R * Math.sin(T) * Math.cos(P);
        VECT[1] = R * Math.sin(T) * Math.sin(P);
        VECT[2] = R * Math.cos(T);
        return VECT;
    }


    public static double eval(final String str) {
        class Parser {
            int pos = -1, c;

            void eatChar() {
                c = (++pos < str.length()) ? str.charAt(pos) : -1;
            }

            void eatSpace() {
                while (Character.isWhitespace(c)) eatChar();
            }

            double parse() {
                eatChar();
                double v = parseExpression();
                if (c != -1) throw new RuntimeException("Unexpected: " + (char) c);
                return v;
            }

            // Grammar:
            // expression = term | expression `+` term | expression `-` term
            // term = factor | term `*` factor | term `/` factor | term brackets
            // factor = brackets | number | factor `^` factor
            // brackets = `(` expression `)`

            double parseExpression() {
                double v = parseTerm();
                for (; ; ) {
                    eatSpace();
                    if (c == '+') { // addition
                        eatChar();
                        v += parseTerm();
                    } else if (c == '-') { // subtraction
                        eatChar();
                        v -= parseTerm();
                    } else {
                        return v;
                    }
                }
            }

            double parseTerm() {
                double v = parseFactor();
                for (; ; ) {
                    eatSpace();
                    if (c == '/') { // division
                        eatChar();
                        v /= parseFactor();
                    } else if (c == '*' || c == '(') { // multiplication
                        if (c == '*') eatChar();
                        v *= parseFactor();
                    } else {
                        return v;
                    }
                }
            }

            double parseFactor() {
                double v;
                boolean negate = false;
                eatSpace();
                if (c == '+' || c == '-') { // unary plus & minus
                    negate = c == '-';
                    eatChar();
                    eatSpace();
                }
                if (c == '(') { // brackets
                    eatChar();
                    v = parseExpression();
                    if (c == ')') eatChar();
                } else { // numbers
                    StringBuilder sb = new StringBuilder();
                    while ((c >= '0' && c <= '9') || (c == '.') || (c == 'E')) {
                        sb.append((char) c);
                        eatChar();
                    }
                    if (sb.length() == 0) throw new RuntimeException("Unexpected: " + (char) c);
                    v = Double.parseDouble(sb.toString());
                }
                eatSpace();
                if (c == '^') { // exponentiation
                    eatChar();
                    v = Math.pow(v, parseFactor());
                }
                if (negate) v = -v; // unary minus is applied after exponentiation; e.g. -3^2=-9
                return v;
            }
        }
        return new Parser().parse();
    }


    public static double randomizeDouble(String OPT, double RANGE) {
        return randomizeDouble(OPT, RANGE, 1E3);
    }

    public static double randomizeDouble(String OPT, double RANGE, double stepShiftInt) {
        Random randomVal = new Random();
        if (OPT.contains("SYM")) {
            return (double) ((randomVal.nextInt((int) stepShiftInt) + 1) / stepShiftInt - .5) * RANGE;
        } else if (OPT.contains("ASYM")) {
            return (double) ((randomVal.nextInt((int) stepShiftInt) + 1) / stepShiftInt) * RANGE;
        }
        return 0;
    }


    public static void randomizeParametersInitial(UnitCell CELL, FragmentData FRAG, int NUM) {

        double dX = 0, dY = 0, dZ = 0;
        double dPhi1 = 0, dPhi2 = 0, dTheta = 0;
        double dO = 0, dU = 0;

        for (int randNUM = 0; randNUM < NUM; randNUM++) {
            for (Iterator<Fragment> iterFRAG = FRAG.getFragMass().iterator(); iterFRAG.hasNext(); ) {
                Fragment itemFRAG = iterFRAG.next();
                dX = dY = dZ = 0;
                dPhi1 = dPhi2 = dTheta = 0;
                dO = dU = 0;
                if (itemFRAG.getFragOPT().contains("X")) {
                    dX = randomizeDouble("SYM", .5);
                }
                if (itemFRAG.getFragOPT().contains("Y")) {
                    dY = randomizeDouble("SYM", .5);
                }
                if (itemFRAG.getFragOPT().contains("Z")) {
                    dZ = randomizeDouble("SYM", .5);
                }


                if (itemFRAG.getFragOPT().contains("R")) {
                    //dPhi1 = randomizeDouble("ASYM", 2 * Math.PI);
                    //dPhi2 = randomizeDouble("ASYM", 2 * Math.PI);
                    //dTheta = randomizeDouble("ASYM", Math.PI);
                    dPhi1 = randomizeDouble("SYM", .5);
                    dPhi2 = randomizeDouble("SYM", .5);
                    dTheta = randomizeDouble("SYM", .5);
                }
                itemFRAG.fragModifyParameters(CELL, dX, dY, dZ, dPhi1, dPhi2, dTheta, dO, dU);
            }
        }
    }

    public static int findNumParameters(FragmentData FRAG) {
        int parametersNum = 0;
        for (Iterator<Fragment> iterFRAG = FRAG.getFragMass().iterator(); iterFRAG.hasNext(); ) {
            Fragment itemFRAG = iterFRAG.next();
            if (itemFRAG.getFragOPT().contains("X")) parametersNum += 1;
            if (itemFRAG.getFragOPT().contains("Y")) parametersNum += 1;
            if (itemFRAG.getFragOPT().contains("Z")) parametersNum += 1;
            if (itemFRAG.getFragOPT().contains("R")) parametersNum += 3;
            if (itemFRAG.getFragOPT().contains("O")) parametersNum += 1;
            if (itemFRAG.getFragOPT().contains("U")) parametersNum += 1;
        }
        return parametersNum;
    }


    public static double[] getParameters(FragmentData FRAG, int[] ParametersList) {
        double[] getPAR = new double[findNumParameters(FRAG)];
        for (int i = 0; i < ParametersList.length; i++) getParameters(FRAG, ParametersList, i);
        return getPAR;
    }


    public static double getParameters(FragmentData FRAG, int[] ParametersList, int numPAR) {
        Fragment itemFRAG = FRAG.getFragMass().get(ParametersList[numPAR] / 10);
        switch (ParametersList[numPAR] % 10) {
            case 1:
                return itemFRAG.getFragX();
            case 2:
                return itemFRAG.getFragY();
            case 3:
                return itemFRAG.getFragZ();
            case 4:
                return itemFRAG.getFragPhi1();
            case 5:
                return itemFRAG.getFragPhi2();
            case 6:
                return itemFRAG.getFragTheta();
            case 7:
                return itemFRAG.getFragO();
            case 8:
                return itemFRAG.getFragU();
        }
        return 0;
    }

    public static void setParameters(FragmentData FRAG, int[] ParametersList, int numPAR, double PAR) {
        Fragment itemFRAG = FRAG.getFragMass().get(ParametersList[numPAR] / 10);
        switch (ParametersList[numPAR] % 10) {
            case 1:
                itemFRAG.setFragX(PAR);
            case 2:
                itemFRAG.setFragY(PAR);
            case 3:
                itemFRAG.setFragZ(PAR);
            case 4:
                itemFRAG.setFragPhi1(PAR);
            case 5:
                itemFRAG.setFragPhi2(PAR);
            case 6:
                itemFRAG.setFragTheta(PAR);
            case 7:
                itemFRAG.setFragO(PAR);
            case 8:
                itemFRAG.setFragU(PAR);
        }
    }

    public static int[] generateParametersList(FragmentData FRAG) {
        int[] gPL = new int[findNumParameters(FRAG)];
        int i = 0, k = 0;
        for (Iterator<Fragment> iterFRAG = FRAG.getFragMass().iterator(); iterFRAG.hasNext(); ) {
            Fragment itemFRAG = iterFRAG.next();
            if (itemFRAG.getFragOPT().contains("X")) {
                gPL[k++] = i * 10 + 1;
            }
            if (itemFRAG.getFragOPT().contains("Y")) {
                gPL[k++] = i * 10 + 2;
            }
            if (itemFRAG.getFragOPT().contains("Z")) {
                gPL[k++] = i * 10 + 3;
            }
            if (itemFRAG.getFragOPT().contains("R")) {
                gPL[k++] = i * 10 + 4;
                gPL[k++] = i * 10 + 5;
                gPL[k++] = i * 10 + 6;
            }
            if (itemFRAG.getFragOPT().contains("O")) {
                gPL[k++] = i * 10 + 7;
            }
            if (itemFRAG.getFragOPT().contains("U")) {
                gPL[k++] = i * 10 + 8;
            }
            i++;
        }
        return gPL;
    }


    public static void addToParameters(UnitCell CELL, FragmentData FRAG, int[] ParametersList, double dPar) {
        for (int i = 0; i < ParametersList.length; i++) addToParameters(CELL, FRAG, ParametersList, i, dPar);
    }

    public static void addToParameters(UnitCell CELL, FragmentData FRAG, int[] ParametersList, double[] dPar) {
        for (int i = 0; i < ParametersList.length; i++)
            addToParameters(CELL, FRAG, ParametersList, i, dPar[i]);
    }

    public static void addToParameters(UnitCell CELL, FragmentData FRAG, int[] ParametersList, int PAR, double dPar) {
        double dX = 0;
        double dY = 0;
        double dZ = 0;
        double dPhi1 = 0;
        double dPhi2 = 0;
        double dTheta = 0;
        double dO = 0;
        double dU = 0;
        Fragment itemFRAG = FRAG.getFragMass().get(ParametersList[PAR] / 10);
        switch (ParametersList[PAR] % 10) {
            case 1:
                dX = dPar;
                break;
            case 2:
                dY = dPar;
                break;
            case 3:
                dZ = dPar;
                break;
            case 4:
                dPhi1 = dPar;
                break;
            case 5:
                dPhi2 = dPar;
                break;
            case 6:
                dTheta = dPar;
                break;
            case 7:
                dO = dPar;
                break;
            case 8:
                dU = dPar;
                break;
            default:
                break;
        }
        itemFRAG.fragModifyParameters(CELL, dX, dY, dZ, dPhi1, dPhi2, dTheta, dO, dU);
    }

    public static List<Integer> integerToRanks(int value){
        List<Integer> result = new ArrayList<>();
        int value_tmp = value;
        int value_goal = 0;
        int counter = 0;
        while(value_goal != value){
            value_tmp %= 10;
            value_goal += Math.pow(10, counter++) * value_tmp;
            result.add(value_tmp);
        }
        return result;
    }

}

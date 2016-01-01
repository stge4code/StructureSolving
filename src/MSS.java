import CrystTools.*;
import Modules.*;

import java.util.Locale;

/**
 * Created by Developer on 10.09.2015.
 */
public class MSS {
    public static void main(String[] args){
        Locale.setDefault(new Locale("en", "US"));

        String FULLPATH = args[0];
        String NAME = args[1];
        String METHOD = args[2];

        System.out.printf("\n%-50s\n", "MODIFIED STRUCTURE SOLUTION v2.1");


        UnitCell CELL = new UnitCell(
                FULLPATH + NAME + ".INS");
        Symmetry SYM = new Symmetry(
                FULLPATH + NAME + ".INS");
        //SYM.printSymmetry();
        FragmentData FRAG = new FragmentData(
                CELL,
                FULLPATH + NAME + ".FGTY",
                FULLPATH + NAME + ".FGLS",
                FULLPATH + NAME + ".OUT",
                FULLPATH + NAME + ".SCAT",
                FULLPATH + NAME + ".ASPC");
        DiffractionData HKL = new DiffractionData(
                CELL,
                FULLPATH + NAME + ".HKL",
                FULLPATH + NAME + ".HKLS");

        PenaltyFunction PSI = new PenaltyFunction(
                FULLPATH + NAME + ".PNLT");

        Energy E = new Energy(HKL, CELL, SYM, PSI);
        //System.out.print(CELL.calcDistance(0, 0, 2));

        if (METHOD.equals("SA")) {
            SAmodule SAROUTINE = new SAmodule(
                    CELL,
                    SYM,
                    FRAG,
                    HKL,
                    E,
                    FULLPATH + NAME + ".SAS",
                    FULLPATH + NAME + ".ARGM",
                    FULLPATH + NAME + ".TEMP",
                    FULLPATH + NAME + ".DUMP");
            SAROUTINE.run();
        }
        if (METHOD.equals("MoR")) {
            MoRmodule MoRROUTINE = new MoRmodule(
                    CELL,
                    SYM,
                    FRAG,
                    HKL,
                    E,
                    FULLPATH + NAME + ".MRS",
                    FULLPATH + NAME + ".TEMP",
                    FULLPATH + NAME + ".DUMP");
            MoRROUTINE.run();
        }

    }

}


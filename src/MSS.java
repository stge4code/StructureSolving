import CrystTools.*;
import Modules.*;

import java.util.Locale;

/**
 * Created by Developer on 10.09.2015.
 */
public class MSS {
    public static void main(String[] args){
        Locale.setDefault(new Locale("en", "US"));

        String FULLPATH = ((args[0]).charAt(args[0].length() - 1) == '\\') ? args[0] : args[0] + '\\';
        String NAME = args[1];
        String METHOD = args[2];

        System.out.printf("\n%-50s\n", "MODIFIED STRUCTURE SOLUTION v2.5");


        UnitCell CELL = new UnitCell(
                FULLPATH + NAME + ".INS");
        Symmetry SYM = new Symmetry(
                FULLPATH + NAME + ".INS");
        //SYM.printSymmetry();
        FragmentData FRAG = new FragmentData(
                CELL,
                FULLPATH + NAME + ".FGTY",
                FULLPATH + NAME + ".STNG",
                FULLPATH  + "_" +  NAME + ".FG",
                FULLPATH + NAME + ".ASCT",
                FULLPATH + NAME + ".ASPC");
        DiffractionData HKL = new DiffractionData(
                CELL,
                FULLPATH + NAME + ".HKL",
                FULLPATH + NAME + ".STNG",
                FULLPATH + "_" + NAME + ".HKL");
        HKL.printHKL();

        PenaltyFunction PSI = new PenaltyFunction(
                FULLPATH + NAME + ".STNG");
        Energy E = new Energy(
                HKL,
                CELL,
                SYM,
                PSI,
                FULLPATH + NAME + ".STNG",
                FULLPATH  + "_" +  NAME + ".HKLF");

        if (METHOD.equals("SA")) {
            SAmodule SAROUTINE = new SAmodule(
                    CELL,
                    FRAG,
                    HKL,
                    E,
                    FULLPATH + NAME + ".STNG",
                    FULLPATH  + "_" +  NAME + ".RF",
                    FULLPATH  + "_" +  NAME + ".DUMP");
            SAROUTINE.run();
        }
        if (METHOD.equals("MoR")) {
            MoRmodule MoRROUTINE = new MoRmodule(
                    CELL,
                    SYM,
                    FRAG,
                    HKL,
                    E,
                    FULLPATH + NAME + ".STNG",
                    FULLPATH  + "_" +  NAME + ".RF",
                    FULLPATH  + "_" +  NAME + ".DUMP");
            MoRROUTINE.run();
        }

    }

}


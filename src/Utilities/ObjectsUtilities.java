package Utilities;

import CrystTools.Atom;
import CrystTools.Fragment;
import CrystTools.SymmetryItem;
import MathTools.FastMath;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by Administrator on 14.09.2015.
 */
public class ObjectsUtilities {
    public static Object deepClone(Object object) {
        try {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ObjectOutputStream oos = new ObjectOutputStream(baos);
            oos.writeObject(object);
            ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
            ObjectInputStream ois = new ObjectInputStream(bais);
            return ois.readObject();
        }
        catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    public static void DeleteFile(String FileName){
        try{
            File file = new File(FileName);
            file.delete();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    public static void saveObject(Object obj2save, String objFileName) throws IOException {
        FileOutputStream fos = new FileOutputStream(objFileName);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(obj2save);
        oos.flush();
        oos.close();
    }
    public static Object recallObject(String objFileName) throws IOException, ClassNotFoundException {
        FileInputStream fis = new FileInputStream(objFileName);
        ObjectInputStream oin = new ObjectInputStream(fis);
        return oin.readObject();
    }
    public static List<String> readFileToList(String filename) {
        List<String> list = new ArrayList<>();
        try (Stream<String> stream = Files.lines(Paths.get(filename))) {
            list = stream.collect(Collectors.toList());
        } catch (IOException e) {
            e.printStackTrace();
        }
        return list;
    }
    public static String generateAtomNum(String name, int num) {
        String charSet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        if (num < 100) {
            return name
                    + String.valueOf(charSet.charAt(num / 10))
                    + String.valueOf(charSet.charAt(num % 10));
        } else {
            //List<Integer> nTor = FastMath.integerToRanks(num);
            if (name.length() < 2) {
                return name
                        + String.valueOf(num / 100)
                        + String.valueOf(charSet.charAt(num / charSet.length()))
                        + String.valueOf(charSet.charAt(num % charSet.length()));
            } else {
                return name
                        + String.valueOf((num % 100))
                        + String.valueOf(charSet.charAt(num % charSet.length()));
            }
        }
    }
    public static <T> void with(T obj, Consumer<T> c) {
        c.accept(obj);
    }

    public static List<String> getContentFromFile(String fileName){
        List<String> result = new ArrayList<>();
        File file = new File(fileName);
        try {
            BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
            try {
                String s;
                while ((s = in.readLine()) != null) {
                    if (!s.isEmpty()) result.add(s);
                }
            } finally {
                in.close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return result;
    }

    public static List<String> getContentFromFile(String fileName, String marker){
        List<String> result = new ArrayList<>();
        File file = new File(fileName);
        try {
            BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
            try {
                String s;
                boolean condition = false;
                while ((s = in.readLine()) != null) {
                    if(s.equals(marker + "{")) {
                        condition = true;
                        continue;
                    }
                    if(condition && s.contains("}")) condition = false;
                    if(condition) {
                        if (!s.isEmpty()) result.add(s);
                    }
                }
            } finally {
                in.close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return result;
    }

    public static void putContentToFile(String fileName, List<String> output){
        putContentToFile(fileName, output, false);
    }
    public static void putContentToFile(String fileName, List<String> output, boolean append){
        File fileOUT = new File(fileName);
        try {
            if (!fileOUT.exists()) {
                fileOUT.createNewFile();
            }

            PrintWriter out = new PrintWriter(new FileWriter(fileOUT.getAbsoluteFile(), append));

            try {
                for(String s : output)  out.println(s);

            } finally {
                out.close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}

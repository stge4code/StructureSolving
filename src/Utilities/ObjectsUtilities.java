package Utilities;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
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
            return name
                    + String.valueOf(num / 100)
                    + String.valueOf(charSet.charAt(num / charSet.length()))
                    + String.valueOf(charSet.charAt(num % charSet.length()));
        }
    }
    public static <T> void with(T obj, Consumer<T> c) {
        c.accept(obj);
    }
}

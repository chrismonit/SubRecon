package subrecon;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class ArrayPrinter {
    
    public static String toString(double[] array, String delim){
        String[] strings = new String[array.length];
        for (int i = 0; i < array.length; i++) {
            strings[i] = Double.toString(array[i]);
        }
        return String.join(delim, strings);
    }
    
    public static void print(double[] array, String delim){
        System.out.println(toString(array,delim));
    }
 
    public static String toString(int[] array, String delim){
       String[] strings = new String[array.length];
        for (int i = 0; i < array.length; i++) {
            strings[i] = Integer.toString(array[i]);
        }
        return String.join(delim, strings);
    }
    
    public static void print(int[] array, String delim){
        System.out.println(toString(array,delim));
    }
    
    public static String toString(Integer[] array, String delim){
       String[] strings = new String[array.length];
        for (int i = 0; i < array.length; i++) {
            strings[i] = Integer.toString(array[i]);
        }
        return String.join(delim, strings);
    }
    
    public static void print(Integer[] array, String delim){
        System.out.println(toString(array,delim));
    }
    
}
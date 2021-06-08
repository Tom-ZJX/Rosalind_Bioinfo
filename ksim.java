import java.io.*;
import java.lang.Math;
import java.util.*;

public class ksim {

    public static String[] ReadFile(String filename) throws Exception {
        File file = new File(filename);
        BufferedReader br = new BufferedReader(new FileReader(file));
        String[] texts = new String[]{"", "", ""};
        String st;
        int index = 1;
        texts[0] = br.readLine();
        while ((st = br.readLine()) != null) {
            texts[index] += st;
            index++;
        }
        return texts;
    }

    public static boolean CheckEditDistance(String s_reconstruct, String t_reconstruct, int k) {
        if (s_reconstruct.length() != t_reconstruct.length()) {
            return false;
        }
        int dist = 0;
        for (int i = 0; i < s_reconstruct.length(); i++) {
            if (s_reconstruct.charAt(i) != t_reconstruct.charAt(i)) {
                dist++;
            }
        }
        return dist <= k;
    }

    public static boolean contains(List<int[]> lst, int[] arr) {
        for (int k = 0; k < lst.size(); k++) {
            if (Arrays.equals(lst.get(k), arr)) {
                return true;
            }
        }
        return false;
    }
    

    public static List<Integer> GetAllIs(String s, String t, int i, int j, int[][] S) {
        if (j == 0) {
            List<Integer> finalIs = new ArrayList<>();
            finalIs.add(i);
            return finalIs;
        }
        int OneMatch = - (s.charAt(i-1) != t.charAt(j-1) ? 1:0);
        List<Integer> finalIs = new ArrayList<>();
        if (S[i][j] == S[i-1][j-1] + OneMatch) {
            finalIs.addAll(GetAllIs(s, t, i-1, j-1, S));
        }
        if (S[i][j] == S[i-1][j] - 1) {
            finalIs.addAll(GetAllIs(s, t, i-1, j, S));
        } 
        if (S[i][j] == S[i][j-1] - 1) {
            finalIs.addAll(GetAllIs(s, t, i, j-1, S));
        }
        return finalIs;
    }

    // try to convert it to non recursive:
    // https://stackoverflow.com/questions/159590/way-to-go-from-recursion-to-iteration
    
    public static List<Integer> GetAllIsIter(String s, String t, int i, int j, int[][] S, int k) {
        Stack<int[]> stack = new Stack<>();
        stack.push(new int[]{i, j, S[i][j]});
        List<Integer> Is = new ArrayList<>();
        Map<String, List<int[]>> map = new HashMap<>();

        while (!stack.empty()) {
            int[] currIndices = stack.pop();
            int currI = currIndices[0];
            int currJ = currIndices[1];
            int score = currIndices[2];
            String key = Integer.toString(currI) + " " + Integer.toString(currJ) + " " + Integer.toString(score);
            if (currJ == 0) {
                if (!Is.contains(currI)) {
                    Is.add(currI);
                }
            } 
            List<int[]> value = map.get(key);
            if (value != null) {
                for (int[] indices: value) {
                    stack.push(indices);
                }
            } else {
                List<int[]> newValue = new ArrayList();

                if (currJ > 0) {
                    int OneMatch = - (s.charAt(currI-1) != t.charAt(currJ-1) ? 1:0);

                    if (S[currI][currJ] == S[currI-1][currJ-1] + OneMatch) {
                        int[] first = new int[]{currI-1, currJ-1, score};
                        newValue.add(first);
                        stack.push(first);
                    }

                    if (S[currI][currJ] == S[currI][currJ-1] - 1) {
                        int[] third = new int[]{currI, currJ-1, score};
                        newValue.add(third);
                        stack.push(third);
                    } else {

                        if (currJ == 1 & S[currI][currJ] == S[currI][currJ-1] & score > -k) {
                            int[] third = new int[]{currI, currJ-1, score-1};
                            newValue.add(third);
                            stack.push(third);
                        }
                    }
                }

                if (currI > 0) {
                    
                    if (S[currI][currJ] == S[currI-1][currJ] - 1) {
                        int[] second = new int[]{currI-1, currJ, score};
                        newValue.add(second);
                        stack.push(second);
                    } else {

                        if (currJ == 0 & S[currI][currJ] == S[currI-1][currJ] & score > -k) {
                            int[] second = new int[]{currI-1, currJ, score-1};
                            newValue.add(second);
                            stack.push(second);
                        }
                    }
                }

            map.put(key, newValue);

        }
                   
    }
    return Is;

}
    
/*
    public static List<Integer> GetAllIsIter(String s, String t, int i, int j, int[][] S, int k) {
        Stack<int[]> stack = new Stack<>();
        stack.push(new int[]{i, j, S[i][j]});
        List<Integer> Is = new ArrayList<>();

        while (!stack.empty()) {
            int[] currIndices = stack.pop();
            int currI = currIndices[0];
            int currJ = currIndices[1];
            int score = currIndices[2];
            if (currJ == 0) {
                if (!Is.contains(currI)) {
                    Is.add(currI);
                }
            }

                    
                if (currJ > 0) {
                    int OneMatch = - (s.charAt(currI-1) != t.charAt(currJ-1) ? 1:0);
                    if (S[currI][currJ] == S[currI-1][currJ-1] + OneMatch) {
                        stack.push(new int[]{currI-1, currJ-1, score});
                    }
                    if (S[currI][currJ] == S[currI][currJ-1] - 1) {
                        stack.push(new int[]{currI, currJ-1, score});
                    } else {

                        if (S[currI][currJ] == S[currI][currJ-1] & score > -k) {
                            stack.push(new int[]{currI, currJ-1, score-1});
                        }
                    }
                }
                    
                if (currI > 0) {
                    if (S[currI][currJ] == S[currI-1][currJ] - 1) {
                        stack.push(new int[]{currI-1, currJ, score});
                    } else {

                        if (S[currI][currJ] == S[currI-1][currJ] & score > -k) {
                            stack.push(new int[]{currI-1, currJ, score-1});
                        }
                    }
                }

                    
                
                
        }
        return Is;

    }
    */

    public static void MotifAlignmentLinear(String s, String t, int k, boolean WriteToFile) throws Exception{
        /**
         * Here I define s as the gene and t as the motif
            k is the max edit distance 
            we want to return the (start loc, length) of all alignments in s with edit distance <= k
            The score 
         */
        int[][] S = new int[s.length()+1][t.length()+1];
        
        //initialization

        for (int j = 0; j < t.length()+1; j++) {
            S[0][j] = -j;
        }
        
        // iteration
        for (int i = 1; i < s.length()+1; i++) {
            for (int j = 1; j < t.length()+1; j++) {
                int OneMatch = - (s.charAt(i-1) != t.charAt(j-1) ? 1:0);
                int[] candidates = new int[] {S[i-1][j-1] + OneMatch, S[i-1][j] - 1, S[i][j-1] - 1};
                Arrays.sort(candidates);
                S[i][j] = candidates[2];
            }
        }

        /*
        for (int i =0; i < s.length()+1; i++) {
            System.out.print(i+": ");
            for (int j = 0; j < t.length()+1; j++) {
                System.out.print(S[i][j] + " "); 
            }
            System.out.println(" ");
        }
        */
        

        List<Integer> Candidates = new ArrayList<>();

        for (int i = 0; i < s.length()+1; i++) {
            if (S[i][t.length()] >= -k) {
                Candidates.add(i);
            }
        } 

        List<String> SubstringLocs = new ArrayList<>();
        // List<String> StringReconstructs = new ArrayList<>();


        for (int m = 0; m < Candidates.size(); m++) {
            int max_i = Candidates.get(m);
            int i = max_i;
            int j = t.length();

            /*
            List<int[]> currIndices = new ArrayList<>();
            currIndices.add(new int[]{i, j});

            List<int[]> Indices = GetIndices(s, t, currIndices, S);

            for (int l = 0; l < Indices.size(); l++) {
                int final_i = Indices.get(l)[0];
                SubstringLocs.add(Integer.toString(final_i+1) + " " + Integer.toString(max_i-final_i));
            }
            */
            
            /*
            while (j > 0) {
                int OneMatch = - (s.charAt(i-1) != t.charAt(j-1) ? 1:0);
                if (S[i][j] == S[i-1][j-1] + OneMatch) {
                    i-=1;
                    j-=1;
                } else {
                    if (S[i][j] == S[i-1][j] - 1) {
                        i-=1;
                    } else {
                        j-=1;
                    }
                }
            }
            */

            /*
            int currScore = S[max_i][t.length()];
            int diff = currScore - (-k);
            for (int temp_i = max_i - t.length() - diff; temp_i < max_i - t.length(); temp_i++) {
                String s_reconstruct = s.substring(temp_i, max_i);
                String t_reconstruct = "-".repeat(s_reconstruct.length() - t.length()) + t;
                if (CheckEditDistance(s_reconstruct, t_reconstruct, k)) {
                    SubstringLocs.add(Integer.toString(temp_i+1) + " " + Integer.toString(max_i-temp_i));
                }

            }
            // this is to account for matches such as 25 5 and 23 7 in ksim_test
            for (int temp_i = Math.max(max_i - t.length(), 0); temp_i < max_i - 1; temp_i++) {
                String t_reconstruct = t;
                String s_reconstruct = s.substring(temp_i, max_i);
                s_reconstruct = "-".repeat(t_reconstruct.length() - s_reconstruct.length()) + s_reconstruct;
                if (CheckEditDistance(s_reconstruct, t_reconstruct, k)) {
                    SubstringLocs.add(Integer.toString(temp_i+1) + " " + Integer.toString(max_i-temp_i));
                }
            }
            */
            List<Integer> Is = GetAllIsIter(s, t, max_i, t.length(), S, k);
            for (int final_i:Is) {
                String newLoc = Integer.toString(final_i+1) + " " + Integer.toString(max_i-final_i);
                if (!SubstringLocs.contains(newLoc)) {
                    SubstringLocs.add(newLoc);
                }
            }
        }

        /*

        for (int m = 0; m < StringReconstructs.size(); m+=2) {
            String s_reconstruct = StringReconstructs.get(m);
            String t_reconstruct = StringReconstructs.get(m+1);
            if (!CheckEditDistance(s_reconstruct, t_reconstruct, k)) {
                System.out.println("Wrong!");
                return;
            }
        }
        */
        
        if (WriteToFile) {
            File myObj = new File("ksim_rslt.txt");
            FileWriter myWriter = new FileWriter("ksim_rslt.txt");
            for (String ss: SubstringLocs) {
                myWriter.write(ss);
                myWriter.write("\n");
            }
            myWriter.close(); // to prevent some part of the text got flushed: https://stackoverflow.com/questions/13426142/bufferedwriter-not-writing-everything-to-its-output-file
        } else {
            for (String ss: SubstringLocs) {
                System.out.println(ss);
            }
        }

    }

    public static void main(String[] args) throws Exception {
        Map<String, int[]> map = new HashMap();
        map.put("1 2 3", new int[]{4,5,6});
        System.out.println(map.get("1 2 3"));
        String[] dnas = new String[] {"ACGGATCGGCATCGT", "ACGTAG"};
        MotifAlignmentLinear(dnas[0], dnas[1], 2, false); // test case

        String[] texts = ReadFile("rosalind_ksim.txt");
        long time = System.nanoTime();
        MotifAlignmentLinear(texts[2], texts[1], Integer.valueOf(texts[0]), true); 
        System.out.println("Alignment Time: " + (System.nanoTime() - time) / 1e9 + " s"); // benchmark time
    }

}

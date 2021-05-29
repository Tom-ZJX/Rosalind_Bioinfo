import java.io.*;
import java.util.Map;
import java.util.HashMap;
import java.lang.Math;
import java.util.Arrays;

public class laff{

    public static String[] ReadFasta(String filename) throws Exception {
        File file = new File(filename);
        BufferedReader br = new BufferedReader(new FileReader(file));
        String[] proteins = new String[]{"", ""};
        String st;
        int index = -1;
        while ((st = br.readLine()) != null) {
            if (st.charAt(0) == '>') {
                index++;
            } else {
                proteins[index] += st;
            }
        }
        return proteins;
    }

    public static void LocalAlignmentAffine(String s, String t, int[][] ScoreMatrix, Map<Character, Integer> m, int StartPenalty, int ExtensionPenalty, boolean WriteToFile) throws Exception {
        int[] Sn = new int[t.length()+1];
        int[] Ss = new int[t.length()+1];
        int[] St = new int[t.length()+1];

        int[][] traceback  = new int[s.length()+1][t.length()+1];

        // initialize all to 3
        for (int i = 0; i < s.length()+1; i++) {
            for (int j = 0; j < t.length()+1; j++) {
                traceback[i][j] = 3; 
            }
        }

        int best = 0;

        int[] OptIndices = new int[]{0, 0};

        int AaToIndexS;
        int AaToIndexT;

        for (int i = 1; i < s.length()+1; i++) {
            int[] new_Sn = new int[t.length()+1];
            int[] new_Ss = new int[t.length()+1];
            int[] new_St = new int[t.length()+1];
            
            AaToIndexS = m.get(s.charAt(i-1));

            for (int j = 1; j < t.length()+1; j++) {
                AaToIndexT = m.get(t.charAt(j-1));
                int curr_score = ScoreMatrix[AaToIndexS][AaToIndexT];

                new_Ss[j] = Math.max(new_Sn[j-1] - StartPenalty,
                                new_Ss[j-1] - ExtensionPenalty);
                new_St[j] = Math.max(Sn[j] - StartPenalty,
                                St[j] - ExtensionPenalty);
                int[] candidates = new int[]{Sn[j-1] + curr_score, new_Ss[j], new_St[j], 0};
                int[] candidates_copy = new int[]{Sn[j-1] + curr_score, new_Ss[j], new_St[j], 0};
                Arrays.sort(candidates);
                new_Sn[j] = candidates[3];

                for (int k = 0; k < 4; k++) {
                    if (candidates_copy[k] == new_Sn[j]) {
                        traceback[i][j] = k;
                    }
                }

                if (new_Sn[j] > best) {
                    best = new_Sn[j];
                    OptIndices[0] = i;
                    OptIndices[1] = j;
                }
            }
                    
            Sn = new_Sn;
            Ss = new_Ss;
            St = new_St;
        }

        // add traceback
        int max_i = OptIndices[0];
        int max_j = OptIndices[1];
        int i = OptIndices[0];
        int j = OptIndices[1];

        while (traceback[i][j] != 3) {
            switch (traceback[i][j]) {
                case 0:
                i-=1;
                j-=1;
                break; // remember to add break everytime

                case 1:
                j-=1;
                break;

                case 2:
                i-=1;
                break;
            }
        }

        if (WriteToFile) {
            File myObj = new File("laff_rslt.txt");
            FileWriter myWriter = new FileWriter("laff_rslt.txt");
            myWriter.write(String.valueOf(best));
            myWriter.write("\n");
            myWriter.write(s.substring(i, max_i));
            myWriter.write("\n");
            myWriter.write(t.substring(j, max_j));
            myWriter.close(); // to prevent some part of the test got flushed: https://stackoverflow.com/questions/13426142/bufferedwriter-not-writing-everything-to-its-output-file

        } else {
            System.out.println(best);
            System.out.println(s.substring(i, max_i)); // s_reconstruct
            System.out.println(t.substring(j, max_j)); // t_reconstruct
        }
            
    }
    
    public static void main(String[] args) throws Exception {

        // creare BLOSUM62 matrix
        File file = new File("BLOSUM62.txt");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String st = br.readLine();
        String[] ar = st.split(" ");
        Map m = new HashMap();
        int counter = 0;
        for (String s: ar) {
            if (s.length() != 0) {
                m.put(s.charAt(0), counter); // Map char to index in BLOSUM62
                counter++;
            }
        }
        int i = 0;
        int j = 0;
        int[][] BLOSUM62 = new int[20][20];

        while ((st = br.readLine()) != null) {
            ar = st.split(" ");
            int index = 1;
            while (index < ar.length) {
                if (ar[index].length() != 0) {
                    BLOSUM62[i][j] = Integer.parseInt(ar[index]);
                    i++;
                }
                index++;
            }
            i = 0;
            j++;
        }
        
        for (int k = 0; k < 20; k++) {
            for (int l = 0; l < 20; l++){
                System.out.print(BLOSUM62[k][l]);
                System.out.print(" ");
            }
            System.out.println(" ");
        } // Verify BLOSUM62 matrix

        long time = System.nanoTime();
        String[] proteins = new String[] {"PLEASANTLY", "MEANLY"};
        LocalAlignmentAffine(proteins[0], proteins[1], BLOSUM62, m, 11, 1, false); // test case
        proteins = ReadFasta("rosalind_laff.txt");
        LocalAlignmentAffine(proteins[0], proteins[1], BLOSUM62, m, 11, 1, true);
        System.out.println("Alignment Time: " + (System.nanoTime() - time) / 1e9 + " s"); // benchmark time
    }
}

import java.io.*;
import java.util.zip.*;


/**
 * 
 */
public class PlotHW {

    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String OUTPUT_DIR = args[1];
        
        int i, j;
        int row, gt00, gt01, gt11, nIndividuals;
        int nCalls_del, nCalls_inv, nCalls_dup, nCalls_ins, nCalls_all;
        String str, gt;
        BufferedReader br;
        BufferedWriter bw_del, bw_inv, bw_dup, bw_ins, bw_all;
        String[] tokens;
        
        // Computing number of individuals
        nIndividuals=0;
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine();
        while (str!=null) {
            if (str.substring(0,6).equalsIgnoreCase("#CHROM")) {
                tokens=str.split("\t");
                nIndividuals=tokens.length-9;
                break;
            }
            str=br.readLine();
        }
        br.close();
        System.err.println(nIndividuals+" individuals");
        
        // Printing GT matrices
        nCalls_del=0; nCalls_inv=0; nCalls_dup=0; nCalls_ins=0; nCalls_all=0;
        bw_del = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/genotypes_del.txt")); bw_del.write("AA,AB,BB\n");
        bw_inv = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/genotypes_inv.txt")); bw_inv.write("AA,AB,BB\n");
        bw_dup = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/genotypes_dup.txt")); bw_dup.write("AA,AB,BB\n");
        bw_ins = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/genotypes_ins.txt")); bw_ins.write("AA,AB,BB\n");
        bw_all = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/genotypes_all.txt")); bw_all.write("AA,AB,BB\n");
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            gt00=0; gt01=0; gt11=0;
            for (i=0; i<nIndividuals; i++) {
                gt=tokens[9+i];
                if ( gt.indexOf("0/0")>=0 || gt.indexOf("./.")>=0 || gt.indexOf("./0")>=0 || gt.indexOf("0/.")>=0 || 
                     gt.indexOf("0|0")>=0 || gt.indexOf(".|.")>=0 || gt.indexOf(".|0")>=0 || gt.indexOf("0|.")>=0
                   ) gt00++;
                else if ( gt.indexOf("0/1")>=0 || gt.indexOf("1/0")>=0 || gt.indexOf("./1")>=0 || gt.indexOf("1/.")>=0 ||
                          gt.indexOf("0|1")>=0 || gt.indexOf("1|0")>=0 || gt.indexOf(".|1")>=0 || gt.indexOf("1|.")>=0
                        ) gt01++;
                else if (gt.indexOf("1/1")>=0 || gt.indexOf("1|1")>=0) gt11++;
            }
            row=svType2Row(getField(tokens[7],SVTYPE_STR));
            if (row==0) {
                nCalls_del++;
                bw_del.write(gt00+","+gt01+","+gt11+"\n");
            }
            else if (row==1) {
                nCalls_inv++;
                bw_inv.write(gt00+","+gt01+","+gt11+"\n");
            }
            else if (row==2) {
                nCalls_dup++;
                bw_dup.write(gt00+","+gt01+","+gt11+"\n");
            }
            else if (row==3) {
                nCalls_ins++;
                bw_ins.write(gt00+","+gt01+","+gt11+"\n");
            }
            nCalls_all++;
            bw_all.write(gt00+","+gt01+","+gt11+"\n");
            str=br.readLine();
        }
        br.close();
        System.err.println("N. calls: "+nCalls_all);
        System.err.println("DEL: "+nCalls_del);
        System.err.println("INV: "+nCalls_inv);
        System.err.println("DUP: "+nCalls_dup);
        System.err.println("INS: "+nCalls_ins);
        bw_del.close(); bw_inv.close(); bw_dup.close(); bw_ins.close(); bw_all.close();
    }
    
    
	private static final int svType2Row(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return 0;
		else if (type.equalsIgnoreCase(INV_STR)) return 1;
        else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return 2;
        else if ( type.equalsIgnoreCase(INS_STR) ||
                  type.equalsIgnoreCase(INS_ME_STR) ||
                  type.equalsIgnoreCase(INS_NOVEL_STR)
                ) return 3;
		else return -1;
	}
    
    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	private static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR)) {
			while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}    
    
    
	/**
	 * Basic constants
	 */
	public static final char COMMENT = '#';
	public static final String SEPARATOR = ";";
	public static final String PASS_STR = "PASS";
	public static final String PRECISE_STR = "PRECISE";
	public static final String IMPRECISE_STR = "IMPRECISE";
	public static final String END_STR = "END";
	public static final String SVTYPE_STR = "SVTYPE";
	public static final String SVLEN_STR = "SVLEN";
	public static final String CHR2_STR = "CHR2";
	public static final String CT_STR = "CT";
	public static final String CT_325_STR = "'3to5'";
	public static final String CT_523_STR = "'5to3'";
	public static final String CT_525_STR = "'5to5'";
	public static final String CT_323_STR = "'3to3'";
    
	/**
	 * SV types: labels used by callers.
	 */
	public static final String DEL_STR = "DEL";
	public static final String DEL_ME_STR = "DEL:ME";
	public static final String DEL_INV_STR = "DEL/INV";
	public static final String INS_STR = "INS";
	public static final String INS_ME_STR = "INS:ME";
	public static final String INS_NOVEL_STR = "INS:NOVEL";
	public static final String DUP_STR = "DUP";
	public static final String DUP_TANDEM_STR = "DUP:TANDEM";
	public static final String DUP_INT_STR = "DUP:INT";
	public static final String INV_STR = "INV";
	public static final String INV_DUP_STR = "INVDUP";
	public static final String CNV_STR = "CNV";
	public static final String BND_STR = "BND";
	public static final String TRA_STR = "TRA";
    
	/**
     *
	 */
	public static final byte TYPE_INSERTION = 1;
	public static final byte TYPE_DELETION = 2;
	public static final byte TYPE_DEL_INV = 3;
	public static final byte TYPE_INVERSION = 4;
	public static final byte TYPE_INV_DUP = 5;
	public static final byte TYPE_DUPLICATION = 6;
	public static final byte TYPE_CNV = 7;
	public static final byte TYPE_BREAKEND = 8;
	public static final byte TYPE_TRANSLOCATION = 9;
    
	/**
	 * Confidence intervals of positions.
	 *
	 * Remark: some callers report a standard deviation instead of a confidence
	 * interval. Sniffles reports additional interval information in its BEDPE
	 * output (an alternative to VCF), which our programs disregard for 
	 * simplicity. Some callers use CILEN to express a "confidence interval 
	 * around inserted/deleted material between breakends": we interpret CILEN
	 * exactly like CIEND, and we ignore it for insertions (since representing
	 * variable-length insertion strings complicates our code).
	 */
	public static final String CI_SEPARATOR = ",";
	public static final String CIPOS_STR = "CIPOS";
	public static final String CIEND_STR = "CIEND";
	public static final String STD_START1_STR = "STD_quant_start";
	public static final String STD_START2_STR = "STD_POS1";
	public static final String STD_END1_STR = "STD_quant_stop";
	public static final String STD_END2_STR = "STD_POS2";
	public static final String CILEN_STR = "CILEN";
	public static final int CIPOS_STR_LENGTH = CIPOS_STR.length();
	public static final int CIEND_STR_LENGTH = CIEND_STR.length();
	public static final int STD_START1_STR_LENGTH = STD_START1_STR.length();
	public static final int STD_START2_STR_LENGTH = STD_START2_STR.length();
	public static final int STD_END1_STR_LENGTH = STD_END1_STR.length();
	public static final int STD_END2_STR_LENGTH = STD_END2_STR.length();
	public static final int CILEN_STR_LENGTH = CILEN_STR.length();
    
}
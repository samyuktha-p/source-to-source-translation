package cetus.transforms;
import cetus.analysis.AnalysisPass;
import cetus.analysis.DDGraph;
import cetus.analysis.DependenceVector;
import cetus.analysis.LoopTools;
import cetus.hir.*;

/*import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.Iterator;*/
import java.util.*;
//import javafx.util.Pair;

import java.io.PrintWriter;

import cetus.analysis.LoopTools;
import cetus.hir.*;
import cetus.application.ChainTools;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

enum DpuCopy{
    DPU_COPY_FROM,
    DPU_COPY_TO,
    DPU_BROADCAST_FROM,
    DPU_BROADCAST_TO
}

enum ArgType {
    GLOBAL_BUFFER,
    NORMAL,
    LOCAL_BUFFER
}

enum AddRef {
    BEFORE,
    AFTER
}

enum AccessType {
    WRITE,
    READ,
    WRITE_READ
}

/*enum PolyStart {
    T_START(0),
    G_START(3),
    B_START(6);
    
    int value;
    
    PolyStart(int value) {
        this.value = value;
    }
    
    int getValue() {
        return value;
    }
}*/

class KernelRepr {  //kernel representation for a kernel_create
    IDExpression cl_kernel_id;
    Expression kernel_name;
    int num_dim = 0;
    int[] order = null;
    
    
    int partition_dim_grid = 1; // 0
    int n_tasklets = 4; // 4
    int multi_wgs = 0;  //-
    
    int partition_dim_wg = 0;   // 0
    int n_multi_wgs = 2;    // 2

    int manual_local_dimensions[] = {2, 2, 2};

    
    WorkHorseArray wha = null;
    
    static boolean POSSIBLE = true;
    boolean local_manual = false;
 
    
//     Map<IDExpression, Expression> parameter_size = new HashMap<IDExpression, Expression>();
    TranslationUnitKernelInfo t_info = null;
    Map<Integer, Expression> parameter_size = new HashMap<Integer, Expression>();
    
    Map<IDExpression, Expression> parameter_offset = new HashMap<IDExpression, Expression>();
    Map<IDExpression, Integer> arguments_map = new HashMap<IDExpression, Integer>();
    Map<Integer, Pair<ArgType, IDExpression>> argument_list = new HashMap<Integer, Pair<ArgType, IDExpression>>();
    Map<Integer, Integer> argument_buffer_map = new HashMap<Integer, Integer>();
    Map<Integer, Pair<Integer, Expression>> index_stride_patterns = new HashMap<Integer, Pair<Integer, Expression>>();

    Long written_buffer_status = (long)0;
    Long read_buffer_status = (long)0;

    CompoundStatement read_comp_block = null;
    Statement read_comp_block_ref = null;
    CompoundStatement while_dpu_block = null;
    
    
//     Map<Integer, Pair<Integer, IDExpression>> read_argument_list = new HashMap<Integer, Pair<Integer, IDExpression>>();

    void set_local_dimensions(boolean val) {
        this.local_manual = val;

    }
    void set_translation_kernel_info(TranslationUnitKernelInfo t_info) {
        this.t_info = t_info;
    }

    TranslationUnitKernelInfo get_translation_kernel_info() {
        return t_info;
    }

    IDExpression get_mapped_global_buffer(int index) {
        Pair<ArgType, IDExpression> pair_gb = argument_list.get(index);
        if(pair_gb.getFirst() == ArgType.GLOBAL_BUFFER)
            return pair_gb.getSecond();
        return null;
    }
    
    IDExpression get_partition_dim_grid() {
        return HostOneTBtoDPU.getIteratorVariable(PolyStart.B_START, partition_dim_grid, "", "l");
    }

    void set_kernel_id(IDExpression id) {
        cl_kernel_id = id;
    }
    
    void set_kernel_name(Expression name) {
        kernel_name = name;
    }

    void set_write_status(int index) {
        this.written_buffer_status = this.written_buffer_status|(1l<<index);

    }
    void set_read_status(int index) {
        this.read_buffer_status = this.read_buffer_status|(1l<<index);
    }

    void setIndexStridePatterns(int arg_index, int num_off_exprs, Expression final_stride) {
        index_stride_patterns.put(arg_index, new Pair(num_off_exprs, final_stride));
    }


    int getNumIndexPatterns(int arg_index) {
        if(index_stride_patterns.containsKey(arg_index))
            return index_stride_patterns.get(arg_index).getFirst();
        return 0;

    }

    Expression getFinalStride(int arg_index) {
        if(index_stride_patterns.containsKey(arg_index))
            return index_stride_patterns.get(arg_index).getSecond();
        return null;
    }

    Long get_read_status() {
        return this.read_buffer_status;
    }
    
    Long get_pending_read_status() {
        return this.written_buffer_status^this.read_buffer_status;
    }

    void put_argument_list(int arg_index, ArgType par_type, IDExpression kernel_arg, GpuBuffers gpu_buffers) {
        argument_list.put(arg_index, new Pair(par_type, kernel_arg));
        arguments_map.put(kernel_arg, arg_index);
        if(par_type == ArgType.GLOBAL_BUFFER) {
            argument_buffer_map.put(arg_index, gpu_buffers.get_index_global_buffer(kernel_arg));
        }
    }

    int get_argument_index(IDExpression id_exp) {
        if(arguments_map.containsKey(id_exp))
            return arguments_map.get(id_exp);
        return -1;
    }

    int get_argument_buffer_map(int arg_index) {
        return argument_buffer_map.get(arg_index);
    }

    Expression get_kernel_name() {
        return kernel_name;
    }
    
    String get_kernel_name_string() {
        return ((StringLiteral)kernel_name).getValue();
    }
    
    WorkHorseArray getWorkHorseArray() {
        return wha;
    }
    
    void set_num_dim(int dim) {
        num_dim = dim;
        if(num_dim<=partition_dim_grid) {
            partition_dim_grid = 0;
            partition_dim_wg = 0;
        }
        order = new int[dim];
        for(int i=0; i<num_dim; i++)
            order[i] = i;
    }
    
    void set_work_horse_array(WorkHorseArray wha) {
        this.wha = wha;
    }
    void add_parameter_size(Integer param_ind, Expression size) {
        parameter_size.put(param_ind, size);
    }
    
    Expression get_parameter_size(Integer param_ind) {
        if(parameter_size.containsKey(param_ind)) {
            return parameter_size.get(param_ind);
        }
        else
            return new IntegerLiteral(0);
    }
    
    void add_parameter_offset(IDExpression param, Expression offset) {
        parameter_offset.put(param, offset);
    }
    
    Expression get_limit_value(IDExpression prim, Primitive primitive, boolean start) {
            int num_args = argument_list.size();
            int primitive_ind = primitive.is_primitive(prim);
            
//             System.out.print("PRIM: " + prim + " " + primitive_ind + " ");
            if(primitive_ind<0) {
                System.out.println("Can't possible to find offset beforehand: " + prim);
            }
            else if(primitive_ind<64) {
                    if(start)
                        return new ArrayAccess(prim.clone(), new IntegerLiteral(0));
                    return new ArrayAccess(prim.clone(), new IntegerLiteral(1));
//                 return wha.getIteratorLimit(wha.convert_to_polystart(primitive_ind), primitive_ind%3);
            }
            else if(primitive_ind==64) {
                return null;
            }
            else if(primitive_ind<(65+num_args)) {
                return argument_list.get(primitive_ind-65).getSecond().clone();
            }
            else if(primitive_ind<96) {
                System.out.print(" total= " + num_args + " ");
            }
            else {
                IDExpression new_prim_name = new NameID("host_"+prim.getName());
                if(start) 
                    return new ArrayAccess(new_prim_name, new IntegerLiteral(0));
                return new ArrayAccess(new_prim_name, new IntegerLiteral(1));
                // RangeExpression re = primitive.iterator_var.get(prim);
                // if(re != null) {
                //     if(start)
                //         return re.getLB();
                //     return re.getUB();
                // }
            }
            return null;
    }
    
    Expression expand_expression(Expression expr, Primitive primitive, boolean start) {
            List<Traversable> children = expr.getChildren();
            Expression dup_expr = expr.clone();
            
            if(dup_expr instanceof IDExpression) {
                Expression new_expr = get_limit_value((IDExpression)dup_expr, primitive, start);
                if(new_expr!=null) {
                    return new_expr;
                }
            }

            int child_num = 0;
            for(Traversable child : children) {
//                 System.out.print(child+ ": ");
                if(child instanceof IDExpression) {
                    IDExpression child_id = (IDExpression) child;
//                     System.out.print(" child_id + ": ");
                    Expression new_expr = get_limit_value(child_id, primitive, start);
                    if(new_expr!=null) {
//                         System.out.println(new_expr);
                        dup_expr.setChild(child_num, new_expr);
                    }
                }
                else if(!(child instanceof FunctionCall) && !(child instanceof Literal) && (child instanceof Expression)) {
                    
                    Expression sub_expr = (Expression) child;
//                         System.out.println(sub_expr.toString());
                    dup_expr.setChild(child_num, expand_expression(sub_expr, primitive, start));
                }
                
                else if(!(child instanceof Literal)){
                    System.out.println("We are regretting that conversion to DPU is not possible : " + child);
                    
                    POSSIBLE = false;
                }
                
                ++child_num;
            }
          return dup_expr;
    }
    
    Expression change_index_values(Expression expr, IntegerLiteral ind) {
        
            Expression dup_expr = expr.clone();
//             List<Traversable> children = expr.getChildren();
            DFIterator<ArrayAccess> df_iterator = new DFIterator<ArrayAccess>(dup_expr, ArrayAccess.class);
            
            int child_num = 0;
//             for(Traversable child : children) {
            while(df_iterator.hasNext()) {
                Traversable child = df_iterator.next();
                
                if(child instanceof ArrayAccess) {
                    ArrayAccess child_id = (ArrayAccess) child;
                   child_id.setIndex(0, ind.clone());
                }
            }
            return dup_expr;
    }
    
    Statement get_offset_expression_end(Expression cpu_sym_id, int offset_ind, Expression exp_offset_expr, String pre) {
//         String post = "_"+Integer.toString(offset_ind);;
//         if(offset_ind == 0) {
//             post = "";
//         }
        String post = "";
        
        Expression lhs_assign_expr = new NameID(pre+cpu_sym_id.toString()+"_end"+post);
        
        Expression assign_expr_offset = new AssignmentExpression(lhs_assign_expr, AssignmentOperator.NORMAL, change_index_values(exp_offset_expr, new IntegerLiteral(1)));
        return new ExpressionStatement(assign_expr_offset);
    }

    static Expression get_offset_expression_size(Expression cpu_sym_id, Expression ele_size_id, String pre, boolean make_align) {
        // String post = "_"+Integer.toString(offset_ind);
        // if(offset_ind == 0) {
        //     post = "";
        // }
        
        // Expression lhs_assign_expr_size = new NameID(pre+cpu_sym_id.toString()+"_size"+post);
        Expression lhs_assign_expr_end = new NameID(pre+cpu_sym_id.toString()+"_end");
        lhs_assign_expr_end = new ArrayAccess(lhs_assign_expr_end, new NameID("copy_i"));
        Expression lhs_assign_expr_offset = new NameID(pre+cpu_sym_id.toString()+"_offset");
        lhs_assign_expr_offset = new ArrayAccess(lhs_assign_expr_offset, new NameID("copy_i"));
//         Expression lhs_assign_expr_offset = new NameID(cpu_sym_id.toString()+"_offset"+post);
        
        
//             lhs_assign_expr_size = new NameID(cpu_sym_id.toString()+"_size_"+Integer.toString(off_ind));
        Expression rhs = new BinaryExpression(lhs_assign_expr_end,BinaryOperator.SUBTRACT, lhs_assign_expr_offset);
        rhs = OneArith.addOne(rhs);
        if(make_align)
            rhs = make8ByteAligned(rhs, ele_size_id.clone());
        
        return rhs;
    }

    // void constructStartEndOffsetStmts(List<Expression> offset_expressions, TranslationUnitKernelInfo t_info, int type, String pre, boolean start) {}

    static Expression make8ByteAligned(Expression arg1, Expression arg2) {
        return new FunctionCall(new NameID("make8ByteAligned"), arg1, arg2);
    }

    void appendOffsetExpressions(List< Triple<Statement, Integer, AccessType> > offset_list, List<Expression> offset_expressions, AccessType acc_type, Long offset_flag, String offset_id_name, int gpu_ind, Primitive primitive, boolean start) {
        int off_ind=0;
        
        for(Expression offset_expr : offset_expressions) {

            if( (gpu_ind!=-1) && ((offset_flag&(((long)1)<<off_ind))==(long)0) ) {
                // if((gpu_ind!=-1) && start && (acc_type == AccessType.READ))
                    // System.out.println("SKIP : (" + off_ind + ") " +  offset_expr);
                off_ind++;
                continue;
            }
            

            Expression lhs_assign_expr_offset= null;


            lhs_assign_expr_offset = new NameID(offset_id_name);

            
            Expression exp_offset_expr = expand_expression(offset_expr, primitive, start);

            if(start && (gpu_ind!=-1)){
                if(acc_type == AccessType.WRITE)
                    System.out.println(offset_id_name+": " +exp_offset_expr+ " - WRITE");
                if(acc_type == AccessType.READ)
                    System.out.println(offset_id_name+": " +exp_offset_expr+ " - READ");
            }
            
            lhs_assign_expr_offset = new ArrayAccess(lhs_assign_expr_offset, new IntegerLiteral(off_ind));
            AssignmentExpression assign_expr_offset = new AssignmentExpression(lhs_assign_expr_offset, AssignmentOperator.NORMAL, exp_offset_expr);

            offset_list.add(new Triple(new ExpressionStatement(assign_expr_offset), off_ind, acc_type) );
            
            off_ind++;
        }        
    }


    
    
    List<Triple<Statement, Integer, AccessType>> get_offset_expression(int param_ind, Expression cpu_sym_id, CompoundStatement dpu_kernel_call_for_block, String pre, boolean start, boolean write_only) {
        //start - is it start / end expression, write_only - do you want indexes which are used for write only or [both read and write]?
        KernelRepr _kernel = this;

        TranslationUnitKernelInfo t_info = _kernel.get_translation_kernel_info();
        Primitive primitive = t_info.primitive;
        IndexToAccessInfo ind_info = t_info.ind_to_access_info.get(param_ind);

        if(ind_info == null)
            return null;
        List<Expression> offset_expressions = ind_info.access_map;
        // if(start && write_only) {
        //     System.out.println("READ: ");
        //     for(Expression expr : offset_expressions) {
        //         System.out.print(", " + expr);
        //     }
        //     System.out.println();
        // }
        Long offset_wflag = ind_info.wmap;
        Long offset_rflag = ind_info.rmap;

        System.out.println("Flag["+param_ind+"]: "+ Long.toBinaryString(offset_wflag)+" "+Long.toBinaryString(offset_rflag));

        Long offset_rwflag = 0l;
        if(offset_wflag!=null && offset_rflag!=null)
            offset_rwflag = offset_wflag&offset_rflag;
        
        if(offset_expressions == null)
            return null;
            
        // if(start) {
            // IDExpression glob_buffer = _kernel.get_mapped_global_buffer(param_ind);
            // if(glob_buffer!=null) {
                // int gpu_ind = gpu_buffers.get_index_global_buffer(glob_buffer);
            int gpu_ind = _kernel.get_argument_buffer_map(param_ind);
            if(start) {
                System.out.println(param_ind + " -> " + gpu_ind + ":: W" + Long.toBinaryString(offset_wflag&(~offset_rwflag)) + ", R&W" + Long.toBinaryString(offset_rwflag) + ", R" + Long.toBinaryString(offset_rflag&(~offset_rwflag)));
            }
            
            if(gpu_ind != -1) {

                // System.out.println("[ERROR] gpu index not found for a global buffer " + glob_buffer);
                if(start) {
                    if(write_only) {
                        _kernel.set_read_status(gpu_ind);
                    }
                    else {
                        _kernel.set_write_status(gpu_ind);
                    }
                }
            }    
        // }
        
//         System.out.println("offset expressions size = " + offset_expressions.size());
        
        String post = "";

        if(start)
            post = "_offset";
        else
            post = "_end";

        String offset_id_name = pre+cpu_sym_id.toString()+post;

        List< Triple<Statement, Integer, AccessType> > offset_list = new ArrayList< Triple<Statement, Integer, AccessType> >();

        appendOffsetExpressions(offset_list, offset_expressions, AccessType.WRITE, offset_wflag&(~offset_rwflag), offset_id_name, gpu_ind, primitive, start);
        appendOffsetExpressions(offset_list, offset_expressions, AccessType.WRITE_READ, offset_rwflag, offset_id_name, gpu_ind, primitive, start);
        
        if(!write_only) {
            // int size_before = offset_list.size();

            appendOffsetExpressions(offset_list, offset_expressions, AccessType.READ, offset_rflag&(~offset_rwflag), offset_id_name, gpu_ind, primitive, start);

            // if(start) {
            //     System.out.println("READ ones are included" + size_before + " -> " + offset_list.size());
            // }
        }

        return offset_list;
    }
}

class GpuBuffers {
    Map<IDExpression, Expression> buffer_size = new HashMap<IDExpression, Expression>();
    Map<IDExpression, Expression> buffer_cpu_map = new HashMap<IDExpression, Expression>();
    List<IDExpression> buffer_entries = new ArrayList<IDExpression>();
    
    void add_buffer_size(IDExpression buffer_id, Expression size) {
        buffer_size.put(buffer_id, size);
        if(buffer_entries.indexOf(buffer_id)==-1)
            buffer_entries.add(buffer_id);
    }
    
    void add_buffer_map(IDExpression buffer_id, Expression cpu_buffer_id) {
        buffer_cpu_map.put(buffer_id, cpu_buffer_id);
        if(buffer_entries.indexOf(buffer_id)==-1)
            buffer_entries.add(buffer_id);
    }

    int get_index_global_buffer(IDExpression buffer_id) {
        return buffer_entries.indexOf(buffer_id);
    }

    IDExpression get_global_buffer_id(int index) {
        return buffer_entries.get(index);
    }
    
    Expression get_buffer_map(IDExpression buffer_id) {
        if(!(buffer_cpu_map.containsKey(buffer_id)))
            return null;
        return buffer_cpu_map.get(buffer_id);
    }
}


class WorkHorseArray {
    Expression work_dim = null;
    int n_dimensions = 0;
    Expression global_work = null;
    Expression local_work = null;
    
    Expression[] global_limits = {null, null, null};
    Expression[] local_limits = {null, null, null};
    
    Integer[] global_limits_int = {1, 1, 1};
    Integer[] local_limits_int = {1, 1, 1};
    
    WorkHorseArray() {
        // for(int i=0; i<3; i++) {
        //     global_limits[i] = new IntegerLiteral(1);
        //     local_limits[i] = new IntegerLiteral(1);
        // }
    }
    
    void set_work_dim(Expression work_dim) {
        this.work_dim = work_dim;
    }

    void set_ndimension(Expression n_dimensions) {
        if(!(n_dimensions instanceof IntegerLiteral))
            return;

        this.n_dimensions = (int)((IntegerLiteral)n_dimensions).getValue();

        for(int dim=this.n_dimensions; dim<3; dim++) {
            global_limits[dim] = new IntegerLiteral(1);
            local_limits[dim] = new IntegerLiteral(1);  
        }

    }
    
    void set_global_limit(int dim, Expression val) {
        if(val==null) {
            System.out.println("CHECK: GLOBAL_DIM not in array literal form");
            return;
        }
        if(dim<3) {
            global_limits[dim] = val;
            if(val instanceof IntegerLiteral)
                global_limits_int[dim] = (int)((IntegerLiteral)val).getValue();
        }
    }
    
    void set_global_limit(Expression[] val) {
        int len = val.length;
        if(val==null) {
            System.out.println("CHECK: GLOBAL_DIM not in array literal form");
            return;
        }
        for(int i=0; i<Math.min(len, 3); i++) {
            if(val[i]==null) {
                System.out.println("CHECK: GLOBAL_DIM not in array literal form " +i);
                return;
            }
            global_limits[i] = val[i];
            global_limits[i].setParens(true);
            if(val[i] instanceof IntegerLiteral)
                global_limits_int[i] = (int)((IntegerLiteral)val[i]).getValue();
        }
    }
    
    Expression get_global_limit(int dim) {
        if(dim<3) {
            if(global_limits[dim]!=null)
                return global_limits[dim].clone();
            else
                return new ArrayAccess(global_work.clone(), new IntegerLiteral((long)dim));
        }
        return null;
    }
    
    void set_local_limit(int dim, Expression val) {
        if(val==null) {
            System.out.println("CHECK: LOCAL_DIM not in array literal form");
            return;
        }
        if(dim<3)
            local_limits[dim] = val;
    }
    
    void set_local_limit(Expression[] val) {
        int len = val.length;
        if(val==null) {
            System.out.println("CHECK: LOCAL_DIM not in array literal form");
            return;
        }
        for(int i=0; i<Math.min(len, 3); i++) {
            if(val[i]==null) {
                System.out.println("CHECK: LOCAL_DIM not in array literal form " +i);
                return;
            }
            local_limits[i] = val[i];
            local_limits[i].setParens(true);
        }
    }
    
    Expression get_local_limit(int dim) {
        if(dim<3) {
            if(local_limits[dim]!=null)
                return local_limits[dim].clone();
            else
                return new ArrayAccess(local_work.clone(), new IntegerLiteral((long)dim));
        }
        return null;
    }
    
    Expression getIteratorLimit(PolyStart start, int dim) {
//             String[] limit_names = new String[]{"ss_tx_end", "ss_ty_end", "ss_tz_end", "ss_gx_end", "ss_gy_end", "ss_gz_end"};
            
            switch(start) {
                case T_START:
                    return get_local_limit(dim);
                    
                case G_START:
                    return get_global_limit(dim);
                    
                case B_START: {
//                     System.out.println(wha.get_global_limit(dim).toString());
//                     System.out.println(wha.get_local_limit(dim).toString());
                    return new BinaryExpression(get_global_limit(dim), BinaryOperator.DIVIDE,  get_local_limit(dim));
                    
                }
            }

            return null;
            
        }
        
        Expression getIteratorLimitExpr(PolyStart start, int dim) {
//             String[] limit_names = new String[]{"ss_tx_end", "ss_ty_end", "ss_tz_end", "ss_gx_end", "ss_gy_end", "ss_gz_end"};
            if(dim>3)
                return null;
                
            String[] limit_names = new String[]{"txl", "tyl", "tzl", "gxl", "gyl", "gzl", "bxl", "byl", "bzl"};
            
            String limit_name = limit_names[start.getValue()+dim];
            

            return new ArrayAccess(new NameID(limit_name), new IntegerLiteral(1));
            
//             return null;
            
        }

        Expression getIteratorStartExpr(PolyStart start, int dim) {
//             String[] limit_names = new String[]{"ss_tx_end", "ss_ty_end", "ss_tz_end", "ss_gx_end", "ss_gy_end", "ss_gz_end"};
            if(dim>3)
                return null;
                
            String[] limit_names = new String[]{"txl", "tyl", "tzl", "gxl", "gyl", "gzl", "bxl", "byl", "bzl"};
            
            String limit_name = limit_names[start.getValue()+dim];
            

            return new ArrayAccess(new NameID(limit_name), new IntegerLiteral(0));
            
//             return null;
            
        }
        
        PolyStart convert_to_polystart(int dim) {
            if(dim<0)
                return null;
            else if(dim<3)
                return PolyStart.T_START;
            else if(dim<6)
                return PolyStart.G_START;
            else if(dim<9)
                return PolyStart.B_START;
            else
                return null;
        }
    
}

public class HostOneTBtoDPU extends TransformPass
{
        public HostOneTBtoDPU(Program program) {
                super(program);
        }

        public String getPassName() {
                return new String("[HostOneTBtoDPU]");
        }
        
        
        
        String dpu_set_name = "ss_dpu_set";
        static String dpu_name = "ss_dpu";
        
        PrintWriter out = new PrintWriter(System.out);
        List<Object> add_swap_stmts;
        List<Statement> ref_swap_stmts;
        List<Statement> comp_swap_stmts;
        List<Boolean> is_swap_next;
        
        List<Object> add_to_stmts;
        List<Statement> ref_to_stmts;
        List<Statement> comp_add_stmts;
        List<Boolean> is_add_next;
        List<Boolean> is_add_before;
        
        List<Statement> rm_to_stmts;
        List<Statement> comp_rm_stmts;
        
        HashSet<Expression> rm_exprs;
        List<Statement> not_rm_stmt;
        List<Expression> ref_swap_exprs;
        
        HashSet<IDExpression> remove_idexp;
        
        Statement cur_stmt = null;
        CompoundStatement cur_comp_stmt = null;
        
        Map<IDExpression, KernelRepr> kernels;
        Map<IDExpression, Expression> Definitions;
        Map<String, TranslationUnitKernelInfo> translation_unit_info;
        
        IDExpression prev_kernel_id = null;
        CompoundStatement prev_while_dpu_block = null;
        CompoundStatement cur_read_comp_block = null;
        Statement cur_read_comp_block_ref = null;
        
        static GpuBuffers gpu_buffers;
        TranslationUnit tran_unit = null;

        // Expression CL_SUCCESS = new IntegerLiteral(1);
        Expression CL_SUCCESS = new NameID("CL_SUCCESS");
        IDExpression DPU_MRAM_OFFSET = new NameID("dpu_mram_offset");
        int global_kernel_count = 0;
        
        public void start()
        {
            translation_unit_info = OneTBtoDPU.translation_unit_info;
            
            DFIterator<TranslationUnit> tran_unit_iter = new DFIterator<TranslationUnit>(program, TranslationUnit.class);
        
        
            
            Procedure proc = null;
            Object obj = null;
            PrintWriter out = new PrintWriter(System.out);
                
            gpu_buffers = new GpuBuffers();
            kernels = new HashMap<IDExpression, KernelRepr>();
            Definitions = new HashMap<IDExpression, Expression>();
            rm_exprs = new HashSet<Expression>();
            not_rm_stmt = new ArrayList<Statement>();
            
            
            while(tran_unit_iter.hasNext()) {
                Traversable tu = tran_unit_iter.next();
                if(tu instanceof TranslationUnit) {
                    tran_unit = (TranslationUnit)tu;
                    String fileName = tran_unit.getInputFilename();
                    String file_extension = "";
                    int ind_ext = fileName.lastIndexOf('.');
                    if (ind_ext > 0) {
                        file_extension = fileName.substring(ind_ext+1);
                    }
                    
                    System.out.println("Current file name: " + fileName + " " + file_extension);
                    if(file_extension.equals("cl")) {
                        System.out.println("SKIP TranslationUnit - HOST");
                        continue;
                    }
//                     if(!file_extension.equals("cpp"))
//                         continue;
                    
                    DepthFirstIterator iter = new DepthFirstIterator(tran_unit);
                    while(iter.hasNext()) {
                obj = iter.next();

                //only manipulting procedures
                
                
                if(obj instanceof Procedure) {
                    proc = (Procedure) obj;
                    System.out.println("Proc Name: "+proc.getName());

                    if(proc.getName().getName().equals("main")) {
                        String code = "FILE *filePointer;\n"+
                            "//static int kernel_num = 0;\n"+
                            "char file_name[30];\n"+
                            "snprintf(file_name, 30, \"dpu_clock_results_%d.txt\", kernel_num);\n"+
                            "filePointer = fopen(file_name, \"a\");\n"+
                            "if(p_val==NULL) {\n"+
                            "fputs(\"\\n\", filePointer);\n"+
                            "}\n"+
                            "else {\n"+
                            "for(int i=0; i<5; i++) {\n"+
                            "fprintf(filePointer, \"%d \", p_val[i]);\n"+
                            "}\n"+
                            "fputs(\"\\n\", filePointer);\n"+
                            "}\n\n"+
                            "//++kernel_num;\n"+
                            "fclose(filePointer);\n";
                        Declaration note_down_decl = createDecls.createProcedure(Specifier.VOID, "note_down", true, "p_val", "kernel_num", createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, note_down_decl);
                        
                        // code = "\nint rem = 8 - (num&7);"
                        //         +"\nrem = (rem==8) ?num :num+rem;"
                        //         +"\nreturn rem/e_size;";

                        code = "\nint rem = 8 - ((num*e_size)&7);"
                                +"\nrem = (rem==8) ?num :num+(rem/e_size);"
                                +"\n//printf(\"Align: %d %d %d\\n\", num, rem, e_size);"
                                +"\nreturn rem;";
                        

                        Declaration align_decl = createDecls.createProcedure(Specifier.INT, "make8ByteAligned", new ArrayList<String>(Arrays.asList("num", "e_size")), null, null, createCodeAnnotStmt(code));;
                        tran_unit.addDeclarationAfter(proc, align_decl);

                        code = "\treturn (a<b) ?a :b;";
                        Declaration proc_decl = createDecls.createProcedure(Specifier.INT, "ss_min", new ArrayList<String>(Arrays.asList("a", "b")), null, null, createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, proc_decl);
                        
                        code = "\treturn (a<b) ?b :a;";
                        proc_decl = createDecls.createProcedure(Specifier.INT, "ss_max", new ArrayList<String>(Arrays.asList("a", "b")), null, null, createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, proc_decl );

                        code = "if(arr[i][0]>arr[j][2] || arr[j][0]>arr[i][2])\n"+
                                "\treturn 0;\n"+

                                "if(arr[i][1]>arr[j][3] || arr[j][1]>arr[i][3])\n"+
                                "\treturn 0;\n"+
                                "return 1;\n";

                        List<Expression> arr_4_expr = new ArrayList<Expression>(Arrays.asList(null, new IntegerLiteral(4)));
                        List<Specifier> arr_4_spec = new ArrayList<Specifier>(Arrays.asList(new ArraySpecifier(arr_4_expr)));
                        proc_decl = createDecls.createProcedure(new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), arr_4_spec, "ss_isIntersect", new ArrayList<String>(Arrays.asList("arr")), null, new ArrayList<String>(Arrays.asList("i", "j")), createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, proc_decl);

                        code = "arr[i][0] = ss_min(arr[i][0], arr[j][0]);\n"+
                                "arr[i][3] = ss_max(arr[i][3], arr[j][3]);\n"+

                                "arr[i][1] = ss_min(arr[i][1], arr[j][1]);\n"+
                                "arr[i][2] = ss_max(arr[i][2], arr[j][2]);\n" +

                                "flag[i] |= flag[j];\n";
                        proc_decl = createDecls.createProcedure(new ArrayList<Specifier>(Arrays.asList(Specifier.VOID)), arr_4_spec, "combineRect", new ArrayList<String>(Arrays.asList("arr")), new ArrayList<String>(Arrays.asList("flag")), new ArrayList<String>(Arrays.asList("i", "j")), createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, proc_decl);

                        code = "    int is_return = 0;\n"+
                                "while(is_return==0) {\n"+
                                //"printf(\"%d\\n\", is_return);\n"+
                                "is_return = 1;\n"+
                                "for(int i=0; i<size; i++) {\n"+
                                "if(map[i]!=i)\n"+
                                    "continue;\n"+
                                "for(int j=i+1; j<size; j++) {\n"+
                                    "if(map[j]!=j)\n"+
                                        "continue;\n"+
                                    "if(ss_isIntersect(arr, i, j)) {\n"+
                                        "combineRect(arr, flag, i, j);\n"+
                                        "map[j]=i;\n"+
                                        "is_return = 0;\n"+
                                    "}\n"+
                                "}\n"+
                            "}\n"+
                        "}\n"+

                        "for(int i=0; i<size; i++) {\n"+
                            "if(map[i]==i)\n"+
                            "\tcontinue;\n"+
                            "int j=i;\n"+
                            "while(map[j]!=j) {\n"+
                                "j = map[j];\n"+
                            "}\n"+
                            "map[i] = j;\n"+

                            "arr[i][0] -= arr[j][0];\n"+
                            "arr[i][1] -= arr[j][1];\n"+

                        "}\n";
                        proc_decl = createDecls.createProcedure(new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), arr_4_spec, "ss_unify", new ArrayList<String>(Arrays.asList("arr")), new ArrayList<String>(Arrays.asList("map", "flag")), new ArrayList<String>(Arrays.asList("size")), createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, proc_decl);

                        code = "int (*arr)[4] = (int (*)[4])malloc(size*sizeof(int)*4);\n"+

                                "for(int i=0; i<size; i++) {\n"+
                                    "arr[i][0] = start[i]%fstride;\n"+
                                    "arr[i][1] = start[i]/fstride;\n"+

                                    "arr[i][2] = end[i]%fstride;\n"+
                                    "arr[i][3] = end[i]/fstride;\n"+

                                    "if(arr[i][0]>arr[i][2])\n"+
                                        "\treturn NULL;\n"+
                                    "if(arr[i][1]>arr[i][3])\n"+
                                        "\treturn NULL;\n"+
                                "}\n"+

                                "ss_unify(arr, map, flag, size);\n"+

                                "for(int i=0; i<size; i++) {\n"+
                                    "if(map[i]==i) {\n"+
                                        "stride[i] = arr[i][2]-arr[i][0]+1;\n"+
                                        "start[i] = arr[i][1]*fstride + arr[i][0];\n"+
                                        "end[i] = arr[i][3]*fstride + arr[i][2];\n"+
                                        "//if(stride[i]<=0)\n"+
                                        "//\tstride[i]=fstride;\n"+
                                    "}\n"+
                                    "//printf(\"%d-\", i);\n"+
                                    "//for(int j=0; j<4; j++)\n"+
                                        "//printf(\"%d \", arr[i][j]);\n"+
                                    "//printf(\"\\n\");\n"+
                                "}\n"+

                                // "for(int i=0; i<size; i++) {\n"+
                                //     "start[i] = arr[i][0] + stride[map[i]]* arr[i][1];\n"+
                                //     "end[i] = arr[i][2] + stride[map[i]]* arr[i][3];\n"+
                                // "}\n"+
                                "return arr;\n";

                        proc_decl = createDecls.createProcedure(new ArrayList<Specifier>(Arrays.asList(Specifier.INT, PointerSpecifier.UNQUALIFIED)), null, "generateRectangle", null, new ArrayList<String>(Arrays.asList("start", "end", "stride", "map", "flag")), new ArrayList<String>(Arrays.asList("fstride", "size")), createCodeAnnotStmt(code));
                        tran_unit.addDeclarationAfter(proc, proc_decl);

                        // FunctionCall func_call = createDpuCopy(new NameID("dpu_buff"), new NameID("dpu_sym_offset"), new NameID("host_buff"), null, new NameID("size_byte"), ArgType.GLOBAL_BUFFER, null, DpuCopy.DPU_COPY_TO);

                        // IfStatement code_if = Strides.generateGeneralisedColAccCopy(new NameID("element_size"), new NameID("host_buff"), new NameID("dpu_buff"), new NameID("ss_temp"), func_call, new NameID("fstride"), null).getFirst();

                        // proc_decl = createDecls.createProcedure(new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), null, "copyToDPU", null, new ArrayList<String>(Arrays.asList("start", "end", "stride", "map", "flag")), new ArrayList<String>(Arrays.asList("fstride", "num_indexes", "element_size")), code_if);
                        // tran_unit.addDeclarationAfter(proc, proc_decl);


                        tran_unit.addDeclarationFirst(createDecls.createCodeAnnotDecl("#include<dpu.h>\n"
                                                        +"#include<inttypes.h>\n"
                                                        +"#include<math.h>\n"
                                                        +"#include<limits.h>\n"
                                                        +"void note_down(int* p_val, int kernel_num);\n"
                                                        +"int make8ByteAligned(int num, int e_size);\n"
                                                        +"int ss_min(int a, int b);\n"
                                                        +"int ss_max(int a, int b);\n"
                                                        +"int* generateRectangle(long long *start, long long *end, long long *stride, long long *map, long long *flag, int fstride, int size);\n"
                                                        +"int ss_unify(int arr[][4], long long * map, long long* flag, int size);\n"
                                                        +"void combineRect(int arr[][4], long long *flag, int i, int j);\n"
                                                        +"int ss_isIntersect(int arr[][4], int i, int j);\n"));
                        // proc.getBody().addDeclaration(createDecls.createCodeAnnotDecl("note_down(-1)"));
                    }

                    CompoundStatement proc_comp_stmt = proc.getBody();
                    DFIterator<Traversable> p_iter = new DFIterator<Traversable>(proc.getBody());
                    
                    cur_stmt = null;
                    cur_comp_stmt = null;
                
                    rm_to_stmts = new ArrayList<Statement>();
                    comp_rm_stmts = new ArrayList<Statement>();
                    
                    add_to_stmts = new ArrayList();
                    ref_to_stmts = new ArrayList<Statement>();
                    comp_add_stmts = new ArrayList<Statement>();
                    is_add_before = new ArrayList<Boolean>();
//                     is_swap_next = new ArrayList<Boolean>();
                    
                    add_swap_stmts = new ArrayList();
                    ref_swap_stmts = new ArrayList<Statement>();
                    comp_swap_stmts = new ArrayList<Statement>();
                    
                    remove_idexp = new HashSet<IDExpression>();
                    ref_swap_exprs = new ArrayList<Expression>();
//                     is_add_next = new ArrayList<Boolean>();
				
                    
//                     Map<IDExpression, Expression> kernel_buffer_size = new HashMap<IDExpression, Expression>();
                    
                    
                    //Map<Symbol, Expression> Definitions = new HashMap<Symbol, Expression>();
                    Map<IDExpression, Expression> Definitions = new HashMap<IDExpression, Expression>();
                    //List<Statement> group_changed_stmts = new ArrayList<Statement>();
                    
				
                    while(p_iter.hasNext()) {
                        Traversable o = p_iter.next();
                        if(o instanceof Statement) {
                            cur_stmt = (Statement) o;
                            
                            if(cur_stmt.getParent() instanceof CompoundStatement){
                                cur_comp_stmt = (CompoundStatement)cur_stmt.getParent();
                            }
                            
                            if(o instanceof DeclarationStatement) {
                                getDefinitions(o, Definitions);
//                                 System.out.println("Def size (decl) = "+ Definitions.size()); 
                            }
                            else if(o instanceof ExpressionStatement) {
                                
                                /*List<Expression> comma_expr_args = new ArrayList<Expression>();
                                comma_expr_args.add(new IntegerLiteral(1));
                                comma_expr_args.add(new NameID("var"));
                                CommaExpression comma_expr = new CommaExpression(comma_expr_args);
                               //cur_comp_stmt.addStatementBefore(cur_stmt, new ExpressionStatement(comma_expr));
                                
                                
                                CompoundStatement comp_dummy = new CompoundStatement();
                                comp_dummy.addStatement(new NullStatement());
                                comp_dummy.annotateBefore(new CodeAnnotation("FOR_EACH "+comma_expr.toString()));
                                if(cur_stmt!=null && cur_comp_stmt!=null)
                                    cur_comp_stmt.addStatementBefore(cur_stmt, comp_dummy);*/
                                
                                Expression expr = ((ExpressionStatement)o).getExpression();
                                
                                if(expr instanceof AssignmentExpression) {
                                    // System.out.println("Expression: " + expr);
                                    getDefinitions(o, Definitions);
//                                     System.out.println("Def size (assignment) = "+ Definitions.size()); 
                                    // if(((AssignmentExpression)expr).getRHS() instanceof FunctionCall)
                                    //     System.out.println("Expression: " + expr);
                                    expressionTranslation(expr, Definitions);
                                }
                                
                                
                                
                                if(expr instanceof FunctionCall) {
                                    FunctionCall fun_call = (FunctionCall) expr;
                                    functionCallTranslation(fun_call);
                                    
                                    /*IDExpression func_name = (IDExpression)fun_call.getName();
                                    String str_func_name = func_name.toString();
                                    if(str_func_name.equals("clCreateProgramWithSource")){
                                        fun_call.swapWith(new IntegerLiteral(0));
                                    }*/
                                }
                                
                                
                            }
                            else if(o instanceof IfStatement) {
                                Expression cond = ((IfStatement)o).getControlExpression();
                                
                                boolean check=true;
                                
                                for(IDExpression id: remove_idexp) {
                                    if(cond.findExpression(id)!=null){
                                        check = false;
                                        break;
                                    }
                                }
                                
                                if(check == false) {
                                    boolean remove_status = false;
                                    DFIterator df_iter_body = new DFIterator(((IfStatement)o).getThenStatement());

                                    while(df_iter_body.hasNext()) {
                                        Object body_obj = df_iter_body.next();
                                        
                                        if(body_obj instanceof CompoundStatement)
                                            continue;
                                        if(body_obj instanceof Statement) {
                                            remove_status = false;

                                            if(body_obj instanceof ReturnStatement) {
                                                remove_status = true;
                                            }

                                            else if(body_obj instanceof ExpressionStatement) {
                                                Expression expr = ((ExpressionStatement)body_obj).getExpression();
                                                if(expr instanceof FunctionCall) {
                                                    FunctionCall func_call = (FunctionCall)expr;
                                                    if(func_call.getName().toString().equals("printf")) {
                                                        // System.out.println("ADDED TO REMOVE: " + func_call);
                                                        remove_status = true;
                                                    }

                                                }

                                                else if(expr instanceof AssignmentExpression) {
                                                    Expression as_expr = ((AssignmentExpression) expr).getRHS();
                                                    if(as_expr instanceof Literal) {
                                                        remove_status = true;
                                                    }
                                                    else if(as_expr instanceof UnaryExpression && ((UnaryExpression)as_expr).getExpression() instanceof Literal)
                                                        remove_status = true; 
                                                    // else
                                                        // System.out.println("IF statements can't be removed: " + as_expr.getClass().getName());  
                                                }
                                            }

                                            if(!remove_status) {
                                                // System.out.println("If statement Can't be removed: "+body_obj);
                                                break;
                                            }
                                            
                                        }                                       
                                        
                                    }

                                    if(remove_status) {
                                        rem_statement(cur_stmt, cur_comp_stmt);
                                    }
                                }
                                // else if(exprEliminator(cond)) {
                                else {
                                    DFIterator<FunctionCall> df_funCall_it = new DFIterator<FunctionCall>(cond, FunctionCall.class);
                                    int num_functions = 0;
                                    while(df_funCall_it.hasNext()) {
                                        FunctionCall func_call = (FunctionCall)df_funCall_it.next();
                                        functionCallTranslation(func_call);
                                        if(exprEliminator(func_call)) {
                                            swap_expression(func_call, null);
                                            // func_call.swapWith(CL_SUCCESS.clone());
                                        }
                                        ++num_functions;
                                    }
                                    // if(num_functions == 0 && cond.findExpression(CL_SUCCESS)!=null)
                                        // rem_statement(cur_stmt, cur_comp_stmt);
                                    
                                    // cond.swapWith(CL_SUCCESS.clone());
                                    // rem_statement(cur_stmt, cur_comp_stmt);
                                }
                                // else {
                                //     DFIterator<FunctionCall> df_funCall_it = new DFIterator<FunctionCall>(o, FunctionCall.class);

                                //     while(df_funCall_it.hasNext()) {
                                //         FunctionCall func_call = (FunctionCall)df_funCall_it.next();
                                //         if(exprEliminator(func_call)) {
                                //             // func_call.swapWith(CL_SUCCESS.clone());
                                //             // rem_statement(cur_stmt, cur_comp_stmt);
                                //             // break;
                                //         }
                                //     }
                                // }
                                
                                
                            }
                        }
                    }
                    /*out.write("*** BEGIN - Def***\n");
                    for(IDExpression key: Definitions.keySet()) {
                        out.write(key.getName()+"= ");
                        Expression e = Definitions.get(key);
                        if(e==null)
                            out.write("null");
                        else					
                            e.print(out);
                        out.write("\n");
                    }
                    out.write("-- END --\n");
                    out.flush();*/
                    
                    /*for(Object stmt_obj: add_to_stmts) {
                        if(stmt_obj instanceof Statement) {
                            proc_comp_stmt.addStatement((Statement)stmt_obj);
                        }
                    }*/
                    add_statements_IR();
                    swap_expressions_IR();
                    swap_statements_IR();
                    
                    rm_statements_IR();
                    
                    //separateDeclInits(proc);
                    
                    
                }
                out.flush();
                    }
                }
            }
            
            
        }  
        
        boolean separateDeclInits(Procedure proc) {
            DepthFirstIterator df_iterator = new DepthFirstIterator(proc);
            CompoundStatement c_st = null;
            
            while(df_iterator.hasNext()) {
                Object obj = df_iterator.next();
                
                if(obj instanceof Statement) {
                    Statement cur_stmt = (Statement)obj;
                    
                    if(cur_stmt.getParent() instanceof CompoundStatement){
                        c_st = (CompoundStatement)cur_stmt.getParent();
                    }
                }
                
                if(obj instanceof DeclarationStatement) {
                    DeclarationStatement decl_stmt = (DeclarationStatement)obj;
                    Declaration decl= decl_stmt.getDeclaration();
                    
                    if(decl instanceof VariableDeclaration) {
                        VariableDeclaration v_decl = (VariableDeclaration)decl;
                        //System.out.println("DECL " + v_decl);
                        
                        List<IDExpression> decl_ids = v_decl.getDeclaredIDs();
                        int total = v_decl.getNumDeclarators();
                        Iterator<IDExpression> iter_id = decl_ids.iterator();
                        
                        for(int index=0; index<total; index++) {
                            VariableDeclarator v_declr = (VariableDeclarator)v_decl.getDeclarator(index);
                            IDExpression v_id = iter_id.next();
                            
                            Initializer init = v_declr.getInitializer();
                            
                            AssignmentExpression asgn_expr = null;
                                
                            if((init!=null)) {
                                Iterator<Traversable> init_value_iterator = init.getChildren().iterator();
                                int init_values = 0, expr_init_values=0;
                                
                                while(init_value_iterator.hasNext()) {
                            
                                    Traversable tt = init_value_iterator.next();
                                    ++init_values;
                                    
                                    Expression init_expr = null;
                                    if(tt instanceof Expression)
                                        init_expr = (Expression)tt;
                                        
                                    if(init_expr instanceof Literal)
                                        continue;
                                    
                                    ++expr_init_values;
                                    
                                    asgn_expr = new AssignmentExpression(v_id.clone(), AssignmentOperator.NORMAL, init_expr.clone());
                                    
                                    
//                                     System.out.println("SEPARATED: "+ decl_stmt);
                                }
                                if(init_values==1 && expr_init_values==1) {
                                    c_st.addStatementAfter(decl_stmt, new ExpressionStatement(asgn_expr));
                                    v_declr.setInitializer(null);
                                }
                                
                                else if(init_values>1 && expr_init_values>0){
                                    init_value_iterator = init.getChildren().iterator();
                                    int init_value = 0;
                                    
                                    while(init_value_iterator.hasNext()) {
                                
                                        Traversable tt = init_value_iterator.next();
                                        
                                        
                                        Expression init_expr = null;
                                        if(tt instanceof Expression)
                                            init_expr = (Expression)tt;
                                            
                                        if(init_expr instanceof Literal) {
                                            ++init_value;
                                            continue;
                                        }
                                        
                                        asgn_expr = new AssignmentExpression(new ArrayAccess(v_id.clone(), new IntegerLiteral(init_value)), AssignmentOperator.NORMAL, init_expr.clone());
                                        c_st.addStatementAfter(decl_stmt, new ExpressionStatement(asgn_expr));
                                        ++init_value;
                                    }
                                    
//                                     v_declr.setInitializer(null);
                                
                                }
                                
                            }
                        
                        }
                    }
                    
                }
            }
            return true;
            
        }
        boolean exprEliminator(Expression expr) {
    
            if(expr instanceof FunctionCall) {
                String[] excpt_list = {"clGetPlatformIDs", "clCreateContextFromType", "clGetContextInfo", "clCreateCommandQueue", "clCreateProgramWithSource", "clBuildProgram", "clReleaseCommandQueue", "clReleaseContext", 
                 "clReleaseMemObject", "clReleaseEvent", "clReleaseProgram", "clGetEventProfilingInfo", "clFinish", "printSupportedImageFormats"};
                List<String> exception_list = Arrays.asList(excpt_list);
                
                IDExpression func_name = (IDExpression)((FunctionCall)expr).getName();
                String str_func_name = func_name.getName();
                if(exception_list.contains(str_func_name)) {
                    return true;
                }
                else if(str_func_name.startsWith("cl_"))
                    return true;
            }
            return false;
        
        }

        Specifier getSpecifierOfIDExpression(IDExpression asgn_id)
        {
            Declaration decl = asgn_id.findDeclaration();
            if(decl instanceof VariableDeclaration) {
                VariableDeclaration v_decl = (VariableDeclaration)decl;
                List<Specifier> spec = v_decl.getSpecifiers();
                if(spec.size()==1)
                    return spec.get(0);
                System.out.println("Specifier: "+ spec.toString()); 
            }
            return null;
        }
        boolean expressionTranslation(Expression expr, Map<IDExpression, Expression> Definitions)
        {
            
            AssignmentExpression asgn_expr = (AssignmentExpression) expr;
            Expression rhs = asgn_expr.getRHS();
            Expression lhs = asgn_expr.getLHS();

            boolean return_status = true;
            IDExpression asgn_id = null;
            
            if(lhs instanceof IDExpression) {
                asgn_id = (IDExpression)lhs;
                
                if(asgn_id==null) {
                    System.out.println("ERROR: assigned one is not idexpression");
                    return false;
                }

                getSpecifierOfIDExpression(asgn_id);

            }
           
            if(rhs instanceof FunctionCall) {

                boolean remove = false;

                FunctionCall fun_call = (FunctionCall) rhs;
                Expression call_name = fun_call.getName();
                String call_name_string = call_name.toString();
                // System.out.println("~ Found: "+ call_name_string);
                
                if(call_name_string.equals("clGetContextInfo")) {
                    Expression err = fun_call.getArgument(fun_call.getNumArguments()-1);
                    if(err instanceof UnaryExpression) {
                        err = ((UnaryExpression)err).getExpression();
                        Statement size_stmt = new ExpressionStatement(new AssignmentExpression(err.clone(), AssignmentOperator.NORMAL, createDecls.createSizeofExpression(new NameID("cl_device_id")) ));
                        add_to_statement(cur_stmt, size_stmt,cur_comp_stmt,true);
                    }

                    remove = true;
                }
                else if(call_name_string.equals("clCreateProgramWithSource")) {
                    return_status = convertCreateProgramWithSource(fun_call);
                }
                else if(call_name_string.equals("clCreateKernel")) {
                    Expression prog_expr = fun_call.getArgument(1);
                    if(prog_expr instanceof IDExpression || prog_expr instanceof StringLiteral) {
                        KernelRepr _kernel = new KernelRepr();
                        kernels.put(asgn_id, _kernel);
                        _kernel.set_kernel_id(asgn_id);
                        
                        if(prog_expr instanceof StringLiteral)
                            _kernel.set_kernel_name(prog_expr);
                        else  {
                            Expression def_kernel_name = Definitions.get((IDExpression)prog_expr);
                            if(def_kernel_name==null)
                                System.out.println("ERROR: Definition of kernel name is not given");
                            _kernel.set_kernel_name(def_kernel_name.clone());
                        }

                        _kernel.set_translation_kernel_info(translation_unit_info.get(_kernel.get_kernel_name_string()));
                                                
                        // System.out.println("def size = "+ Definitions.size());
                        
                        /*for(Map.Entry<IDExpression, Expression> def: Definitions.entrySet()) {
                            System.out.println("CURRENT DEF: "+ def.getKey() + " " + def.getValue());
                        }*/
                        // System.out.println("Kernel Name - set " + prog_expr + " " + Definitions.get((IDExpression)prog_expr));
                        /*Declaration decl = ((IDExpression)prog_expr).findDeclaration();
                        if(decl instanceof VariableDeclaration) {
                            VariableDeclaration v_decl = (VariableDeclaration)decl;
                            
                            //init_value(proc, prog_expr);
                        }*/
                    }
                    else {
                        System.out.println("Error: clCreateKernel: Wrong parameters : "+prog_expr.getClass().getName() + " "+ prog_expr);
                        return false;
                    }
                    remove = true;
                }
                else if(call_name_string.equals("clCreateBuffer")) {
                    Expression b_size = fun_call.getArgument(2); 
                    gpu_buffers.add_buffer_size(asgn_id, b_size);
                    Expression err = fun_call.getArgument(fun_call.getNumArguments()-1);
                    if(err instanceof UnaryExpression) {
                        err = ((UnaryExpression)err).getExpression();
                        lhs.swapWith(err.clone());
                    }
                    remove = true;
                    
                }
                else if(call_name_string.equals("clEnqueueWriteBuffer")) {          //CAUTION: below two calls present at two places
                    return_status = convertEnqueueWriteBuffer(fun_call);
                    // remove_idexp.add(asgn_id);
                    remove = true;
                
                }
                else if(call_name_string.equals("clEnqueueReadBuffer")) {          //CAUTION: below two calls present at two places
                    System.out.println("expressionTranslation: "+rhs);
                    return_status = convertEnqueueReadBuffer(fun_call);
                    remove = true;
                
                }
                else if(call_name_string.equals("clEnqueueNDRangeKernel")) {
                    return_status = convertEnqueueNDRangeKernel(fun_call);
                    System.out.println("expressionTranslation: "+rhs);
                }
                else if(call_name_string.equals("clSetKernelArg")) {
                    return_status = convertSetKernelArg(fun_call);
                    remove = true;
                }
                else {
                    remove = exprEliminator(rhs);
                }
                
                
                if(remove) {
//                     cur_stmt.swapWith(new NullStatement());
                    // rem_statement(cur_stmt, cur_comp_stmt);
                    System.out.println("expressionTranslation: "+rhs);

                    rhs.swapWith(CL_SUCCESS.clone());
                    remove_idexp.add(asgn_id);
                    
                }
                // else
                    // System.out.println("~ Encountered: "+ call_name_string);
            }
            return return_status;
        }

        boolean convertCreateProgramWithSource(FunctionCall fun_call){
            // System.out.println("Encountered: clCreateProgramWithSource\n");
            Expression err = fun_call.getArgument(fun_call.getNumArguments()-1);
            if(err instanceof UnaryExpression) {
                err = ((UnaryExpression)err).getExpression();
                if(err instanceof IDExpression)
                    remove_idexp.add((IDExpression)err);
            }
            swap_expression(fun_call, null);
            return true;
        }

        boolean convertEnqueueWriteBuffer(FunctionCall fun_call) {          // clEnqueueWriteBuffer
            IDExpression gpu_buf = (IDExpression)fun_call.getArgument(1);
            Expression cpu_buf = (Expression)fun_call.getArgument(5);
            if(gpu_buffers == null) {
                System.out.println("ERROR: convertEnqueueWriteBuffer(): gpu_buffers are not yet created");
                return false;
            }
            gpu_buffers.add_buffer_map(gpu_buf, cpu_buf);
            return true;
        }
        
        boolean convertEnqueueReadBuffer(FunctionCall fun_call) {           // clEnqueueReadBuffer
            IDExpression gpu_buf = (IDExpression)fun_call.getArgument(1);
            Expression cpu_buf = (Expression)fun_call.getArgument(5);

            if(gpu_buffers.get_buffer_map(gpu_buf)==null) {
                gpu_buffers.add_buffer_map(gpu_buf, cpu_buf);
                System.out.println("Buffer Mapping Added: "+ gpu_buf + " "+ cpu_buf);
            }
//                     gpu_buffers.add_buffer_map(gpu_buf, cpu_buf);
            // remove_idexp.add(asgn_id);
            int k_ind = -1;

            for(Map.Entry<IDExpression, KernelRepr> kernel: kernels.entrySet()) {
                ++k_ind;
                KernelRepr _kernel = kernel.getValue(); // kernels.get(prev_kernel_id);

                 if(_kernel == null || _kernel.getWorkHorseArray()==null) {
                    System.out.println("ERROR: convertEnqueueReadBuffer(): Either Kernel or Work Horse Dim not yet created");
                    if(_kernel==null) {
                        System.out.println("[ISSUE] kernel");
                        return false;
                    }
                    else
                        System.out.println("[ISSUE] WorkHorseArray");
                }

                if(gpu_buffers.get_index_global_buffer(gpu_buf)==-1)
                    System.out.println("[ERROR] gpu buffer not found " + gpu_buf);

                System.out.println("\npending read & current buffer: (" + _kernel.get_kernel_name() + ") " + Long.toBinaryString(_kernel.get_pending_read_status())  + " " + Long.toBinaryString((1l<<gpu_buffers.get_index_global_buffer(gpu_buf))));

                if((_kernel.get_pending_read_status()&(1l<<gpu_buffers.get_index_global_buffer(gpu_buf)))==0) {
                    continue;
                }

                cur_read_comp_block = _kernel.read_comp_block;

               
                if(cur_read_comp_block==null) {
                    cur_read_comp_block = _kernel.read_comp_block = new CompoundStatement();

                    List<Statement> one_block_inner_stmts = createOneBlockLimits(_kernel.getWorkHorseArray(), _kernel.order, true, "p_", _kernel.partition_dim_grid);
                    for(Statement st: one_block_inner_stmts) {
                        if(st instanceof DeclarationStatement) {
                            cur_read_comp_block.addDeclaration(((DeclarationStatement)st).getDeclaration().clone());
                        }
                        else
                            cur_read_comp_block.addStatement(st);
                    }
                    
                    List<Statement> inner_stmts = createIncBlockLimits(_kernel.getWorkHorseArray(), _kernel.order, "p_", _kernel.partition_dim_grid);
                    _kernel.read_comp_block_ref = inner_stmts.get(1);
    //                         Statement out_of_bound_stmt = createifLimitCondition("p_");
                    cur_read_comp_block.addStatement(createifLimitCondition("p_"));

                    for(Statement st: inner_stmts) {
                        if(st instanceof DeclarationStatement) {
                            cur_read_comp_block.addDeclaration(((DeclarationStatement)st).getDeclaration().clone());
                        }
                        else
                            cur_read_comp_block.addStatement(st);
                    }
                    
                    _kernel.while_dpu_block.addStatement(createDpuForEach());
                    _kernel.while_dpu_block.addStatement(cur_read_comp_block);
                    _kernel.while_dpu_block.addStatement(createCountCyclesDpu(k_ind));
                }
                
                convertReadBuffer(fun_call, _kernel);
            }
//             if(cur_read_comp_block == null) {
//                 cur_read_comp_block = new CompoundStatement();

//                 KernelRepr _kernel = kernels.get(prev_kernel_id);
//                 cur_read_comp_block = _kernel.read_comp_block;

//                 if(_kernel == null || _kernel.getWorkHorseArray()==null) {
//                     System.out.println("ERROR: convertEnqueueReadBuffer(): Either Kernel or Work Horse Dim not yet created");
//                     return false;
//                 }
                
//                 List<Statement> one_block_inner_stmts = createOneBlockLimits(_kernel.getWorkHorseArray(), _kernel.order, true, "p_");
//                 for(Statement st: one_block_inner_stmts) {
//                     if(st instanceof DeclarationStatement) {
//                         cur_read_comp_block.addDeclaration(((DeclarationStatement)st).getDeclaration().clone());
//                     }
//                     else
//                         cur_read_comp_block.addStatement(st);
//                 }
                
//                 List<Statement> inner_stmts = createIncBlockLimits(_kernel.getWorkHorseArray(), _kernel.order, "p_");
//                 cur_read_comp_block_ref = inner_stmts.get(0);
// //                         Statement out_of_bound_stmt = createifLimitCondition("p_");
//                 cur_read_comp_block.addStatement(createifLimitCondition("p_"));

//                 for(Statement st: inner_stmts) {
//                     cur_read_comp_block.addStatement(st);
//                 }
                
//                 _kernel.while_dpu_block.addStatement(createDpuForEach());
//                 _kernel.while_dpu_block.addStatement(cur_read_comp_block);
//                 _kernel.while_dpu_block.addStatement(createCountCyclesDpu());
//             }
            
            // convertReadBuffer(fun_call, kernels.get(prev_kernel_id));
            
//                     cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, convertReadBuffer(fun_call, kernels.get(prev_kernel_id)));
//                     cur_read_comp_block.addStatement(convertReadBuffer(fun_call));
            
//                     swap_with_statement(cur_stmt, , cur_comp_stmt);
            return true;
        }
        
        boolean convertEnqueueNDRangeKernel(FunctionCall fun_call) {       // clEnqueueNDRangeKernel    CAUTION: check if it's the only statement in that line (not in declaraiotion stmt with multilple declarators)
            IDExpression kernel_id = (IDExpression)fun_call.getArgument(1);
            KernelRepr _kernel = kernels.get(kernel_id);
            
            if(_kernel==null) {
                System.out.println("ERROR: convertEnqueueNDRangeKernel(): kernel not yet created");
                return false;
            }

            CompoundStatement dpu_block = createDpuKernel(fun_call, _kernel);
            CompoundStatement dpu_block_init = createDpuKernelInit(_kernel);
//                     int[] order = {0, 1, 2};
            int[] order = kernels.get(kernel_id).order;
            WhileLoop while_dpu_block = createWhileBlockLoop(_kernel.wha, order, dpu_block, true, _kernel.partition_dim_grid);
            Statement whole_dpu_block = createBlock(while_dpu_block);
            // CompoundStatement whole_dpu_block = _kernel.whole_dpu_block;
            // whole_dpu_block.addStatement(while_dpu_block);

            _kernel.while_dpu_block = (CompoundStatement)while_dpu_block.getBody();

            //add_to_stmts.add(dpu_block);
            
            CompoundStatement ndrange_block = new CompoundStatement();
            
            ndrange_block.addStatement(createAllocDpu());
            ndrange_block.addStatement(createDecls.createCodeAnnotStmt("\n"));

            
            ndrange_block.addStatement(createDpuLoad(_kernel.get_kernel_name()));
            ndrange_block.addStatement(createDecls.createCodeAnnotStmt("\n"));
            createDpuParamDeclarations(_kernel, ndrange_block, order);  // offset, size decls
            ndrange_block.addStatement(createDecls.createCodeAnnotStmt("\n"));

            dpu_block_init.addStatement(makeDpuAssert(createDpuCopy(new NameID("GRID_LIMIT"), new IntegerLiteral(0), getIteratorVariable(PolyStart.B_START, _kernel.partition_dim_grid, "limit"), null, null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO)));

            createDpuParamCopy(_kernel, ndrange_block, dpu_block_init);      //attaches the parameter variables dpu copies like ntasklets, number of blocks to each dpu...
            
            // ndrange_block.addStatement(createDpuParamTrace());

            ndrange_block.addStatement(createDpuForEach());
            ndrange_block.addStatement(dpu_block_init);
            ndrange_block.addStatement(createDpuLaunch());
            ndrange_block.addStatement(createDecls.createCodeAnnotStmt("\n"));
            
            ndrange_block.addStatement(whole_dpu_block);
            
            ndrange_block.addStatement(createDecls.createCodeAnnotStmt("\n"));
            //ndrange_block.addStatement(createCountCyclesDpu());
            ndrange_block.addStatement(createFreeDpu());
            
            swap_with_statement(cur_stmt, ndrange_block, cur_comp_stmt);
                    
            // prev_kernel_id = kernel_id;         //for readbuffer statement
            // cur_read_comp_block = null;
            
            // not_rm_stmt.add(cur_stmt);
            
            // prev_whole_dpu_block = (CompoundStatement) whole_dpu_block;
            return true;
        }

        boolean convertSetKernelArg(FunctionCall fun_call) {
            Expression kernel_id = fun_call.getArgument(0);
            Expression param_num = fun_call.getArgument(1);
            Expression param_id = fun_call.getArgument(3);
            Expression param_size = fun_call.getArgument(2);
            
            Integer arg_index = (int)((IntegerLiteral)param_num).getValue();
            KernelRepr _kernel = kernels.get((IDExpression)kernel_id);
            if(_kernel==null) {
                System.out.println("ERROR: functionCallTranslation: convertSetKernelArg(): Kernel is not created");
                return false;
            }
            _kernel.add_parameter_size(arg_index, param_size);
            
            DFIterator<Traversable> expr_dfs = new DFIterator<Traversable>(param_id);
            boolean id_present= false;
            while(expr_dfs.hasNext()) {
                //System.out.println(": "+expr_dfs.next().getClass().getName());
                Traversable tt = expr_dfs.next();
                if(tt instanceof IDExpression) {
                    id_present = true;
                    IDExpression kernel_arg = (IDExpression)tt;
                    if(gpu_buffers.buffer_size.containsKey(kernel_arg)) {
                        // System.out.println(arg_index + ": Kernel Global Buffer");    - 0
                        _kernel.put_argument_list(arg_index, ArgType.GLOBAL_BUFFER, kernel_arg, gpu_buffers);
                    }
                    else {
                        // System.out.println(arg_index + " Scalar");   - 1
                        _kernel.put_argument_list(arg_index, ArgType.NORMAL, kernel_arg, gpu_buffers);
                    }
                }    
            }


            if(!id_present) {
                // System.out.println(arg_index + ": Kernel Local buffer");     -2
                _kernel.put_argument_list(arg_index, ArgType.LOCAL_BUFFER, null, gpu_buffers);
            }
            // System.out.println("-------");
            return true;
        }
        
        boolean functionCallTranslation(FunctionCall fun_call) {
            
            boolean remove = false, swap = false, return_status =  true;
        
            Expression call_name = fun_call.getName();
            String call_name_string = call_name.toString();
            if(call_name_string.equals("clSetKernelArg")) {
                return_status = convertSetKernelArg(fun_call);
                remove = true;
            }

            else if(call_name_string.equals("clCreateProgramWithSource")) {
                return_status = convertCreateProgramWithSource(fun_call);
                System.out.println("functionCallTranslation(): "+fun_call);
            }

            else if(call_name_string.equals("clEnqueueWriteBuffer")) {          //CAUTION: below two calls present at two places
                return_status = convertEnqueueWriteBuffer(fun_call);
                swap = true;

                // remove = true;
                
            }
            else if(call_name_string.equals("clEnqueueReadBuffer")) {          //CAUTION: below two calls present at two places
                return_status = convertEnqueueReadBuffer(fun_call);
                swap = true;

                // remove = true;
            }   
            else if(call_name_string.equals("clEnqueueNDRangeKernel")) {
                return_status = convertEnqueueNDRangeKernel(fun_call);
                System.out.println("functionCallTranslation(): "+fun_call);
            }
            else {
                swap = exprEliminator(fun_call);
            }

            // if(2==2)
                // return false;
            if(remove) {
//                 cur_stmt.swapWith(new NullStatement());
                // System.out.println("PARENT: " + fun_call.getParent().getClass().getName() + " " + (fun_call.getParent() == cur_comp_stmt));
                // if(fun_call.getParent() == cur_comp_stmt)
                System.out.println("functionCallTranslation(): "+fun_call);
                rem_statement(cur_stmt, cur_comp_stmt);
                // else
                    // fun_call.swapWith(CL_SUCCESS);                                                   //CHECK: directly executing this is changing AST and producing duplicate results
            }
            else if(swap) {
                System.out.println("functionCallTranslation(): "+fun_call);
                swap_expression(fun_call, null);
            }
            return true;
        }
        
        Expression init_value(Traversable trav, Expression expr) {
            if(expr instanceof IDExpression) {
                Declaration decl = ((IDExpression)expr).findDeclaration();
                if(decl instanceof VariableDeclaration) {
                    VariableDeclaration v_decl = (VariableDeclaration)decl;
                    //List<IDExpression> decl_ids = ;
                    int index = v_decl.getDeclaredIDs().indexOf((IDExpression)expr);
                    if(index!=-1) {
                        VariableDeclarator v_declr = (VariableDeclarator)v_decl.getDeclarator(index);
                        Initializer init = v_declr.getInitializer();
                        if(init!=null) {
                            //out.write(id_sym.getSymbolName() + " " + init.getClass().getName());
                            Iterator<Traversable> init_value_iterator = init.getChildren().iterator();
                            while(init_value_iterator.hasNext()) {
                                
                                Traversable tt = init_value_iterator.next();
                                //out.write(tt.getClass().getName()+" ");
                                if(tt instanceof Expression) {
                                    Statement stmt = ((Expression) tt).getStatement();
                                    //System.out.println("Statement:" + stmt.toString());
                                    return (Expression) tt;
                                }
                            }
                        }
                    }
                }
            }
                
            Expression first_expr = IRTools.findExpression(trav, expr);
            Statement stmt = expr.getStatement();
            // System.out.println("Statement:" + stmt.toString());
            return first_expr;
            
        }
        
        void getDefinitions(Traversable trav, Map<IDExpression, Expression> Definitions) {
            
            if(trav instanceof DeclarationStatement) {
                Declaration decl = ((DeclarationStatement) trav).getDeclaration();
                if(decl instanceof VariableDeclaration) {
                    VariableDeclaration v_decl = (VariableDeclaration)decl;
                    //List<Symbol> decl_syms = v_decl.getDeclaredSymbols();
                    List<IDExpression> decl_ids = v_decl.getDeclaredIDs();
                    int total = v_decl.getNumDeclarators();
                    Iterator<IDExpression> iter_id = decl_ids.iterator();
                    for(int index=0; index<total; index++) {
                        VariableDeclarator v_declr = (VariableDeclarator)v_decl.getDeclarator(index);
                        IDExpression v_id = iter_id.next();
                        Initializer init = v_declr.getInitializer();
                        if(init!=null) {
                            //out.write(id_sym.getSymbolName() + " " + init.getClass().getName());
                            Iterator<Traversable> init_value_iterator = init.getChildren().iterator();
                            int multi_value_counter = 0;
                            while(init_value_iterator.hasNext()) {
                                
                                Traversable tt = init_value_iterator.next();
//                                 System.out.println(v_id.getName()+" " );
                                //if(v_id.getName().equals("global_work"))
                                //    System.out.println("Global_work: "+tt.getClass().getName()+" ");
                                if(tt instanceof Expression) {
                                    Statement stmt = ((Expression) tt).getStatement();
//                                     System.out.println((Expression)tt);
//                                     System.out.println("Statement:" + stmt.toString());
                                    Definitions.put(v_id, (Expression)tt);
                                    ++multi_value_counter;
                                    // break;

                                }
                            }

                            if(multi_value_counter>1) {
                                Definitions.put(v_id, null);
                            }
                            else if(multi_value_counter==1) {
                                Expression tt = (Expression)init.getChildren().iterator().next();
                                if(tt instanceof FunctionCall)
                                    functionCallTranslation((FunctionCall)tt);
                            }

                        }
                        else {
                            Definitions.put(v_id, null);
                        }
                    }
                }
            }
            else if(trav instanceof ExpressionStatement && ((ExpressionStatement)trav).getExpression() instanceof AssignmentExpression) {
                AssignmentExpression ass_expr = (AssignmentExpression)((ExpressionStatement)trav).getExpression();
                Expression lhs = ass_expr.getLHS();	
                if(lhs instanceof IDExpression) {
                    //Set<Symbol> lhs_syms = SymbolTools.getAccessedSymbols(lhs);
                    IDExpression v_id = (IDExpression)lhs;				
                    Definitions.put(v_id, ass_expr.getRHS());
                    //int local_dep = getDependenceVector(ass_expr.getRHS(), dep_vector, proc);
                    //dep_vector.put(lhs_sym, local_dep);
                }
            }
        }
        
        Expression[] getExpansion(IDExpression id_exp) {
            Declaration decl = id_exp.findDeclaration();
            Expression[] ret_expr = new Expression[3];
            if(decl instanceof VariableDeclaration) {
                VariableDeclaration v_decl = (VariableDeclaration)decl;
                //List<Symbol> decl_syms = v_decl.getDeclaredSymbols();
                List<IDExpression> decl_ids = v_decl.getDeclaredIDs();
                int total = v_decl.getNumDeclarators();
                Iterator<IDExpression> iter_id = decl_ids.iterator();
                
                for(int index=0; index<total; index++) {
                    IDExpression v_id = iter_id.next();
                    if(!v_id.equals(id_exp))
                        continue;
                    
                    VariableDeclarator v_declr = (VariableDeclarator)v_decl.getDeclarator(index);
                    
                    Initializer init = v_declr.getInitializer();
                    if(init!=null) {
                        //out.write(id_sym.getSymbolName() + " " + init.getClass().getName());
                        Iterator<Traversable> init_value_iterator = init.getChildren().iterator();
                        int ind = 0;
                        while(init_value_iterator.hasNext()) {
                            
                            Traversable tt = init_value_iterator.next();
                            //System.out.println(v_id.getName()+" ");
//                             if(v_id.getName().equals("global_work"))
//                                 System.out.println("Global_work: "+tt.getClass().getName()+" ");
                            if(tt instanceof Expression) {
//                                 Statement stmt = ((Expression) tt).getStatement();
    //                                     System.out.println("Statement:" + stmt.toString());
//                                 Definitions.put(v_id, (Expression)tt);
                                ret_expr[ind++] = (Expression)tt;
                                //return (Expression) tt;
                            }
                        }
                    }
                }
            }
            return ret_expr;
        }
        
        static Statement makeDpuAssert(Expression expr) {
            IDExpression func_name = new NameID("DPU_ASSERT");
            return new ExpressionStatement( new FunctionCall(func_name, expr) );
        }


        Statement createElementSizeStmt(Expression cpu_buffer_id, Expression ele_size_id)
        {
            Expression ele_size_expr = null;
            if(SymbolTools.isPointer(cpu_buffer_id)) {
                List<Specifier> specs = new ArrayList(Arrays.asList(Specifier.CHAR, PointerSpecifier.UNQUALIFIED));
                Expression _lhs = new BinaryExpression(cpu_buffer_id.clone(), BinaryOperator.ADD, new IntegerLiteral(1));
                Expression _rhs = cpu_buffer_id.clone();
                ele_size_expr = new BinaryExpression(new Typecast(specs, _lhs), BinaryOperator.SUBTRACT, new Typecast(specs, _rhs));
            }
            else
                ele_size_expr = new BinaryExpression(createDecls.createSizeofExpression(cpu_buffer_id.clone()), BinaryOperator.DIVIDE, createDecls.createSizeofExpression(new ArrayAccess(cpu_buffer_id.clone(), new IntegerLiteral(0)) ));

            ExpressionStatement ele_size = new ExpressionStatement(new AssignmentExpression(ele_size_id.clone(), AssignmentOperator.NORMAL, ele_size_expr));
            return ele_size;           
        }

        Statement createElementSizeStmt(Expression ele_size_id, Specifier spec) {
            Expression ele_size_expr = null;
            ele_size_expr = new SizeofExpression(new ArrayList<Specifier>(Arrays.asList(spec)) );
            ExpressionStatement ele_size = new ExpressionStatement(new AssignmentExpression(ele_size_id.clone(), AssignmentOperator.NORMAL, ele_size_expr));
            return ele_size;
        }
        
        CompoundStatement createDpuKernel(FunctionCall kernelND, KernelRepr _kernel) {  //CAUTION: global_work, local_work should be initiazed in the declaration statement
            if(_kernel==null) {
                System.out.println("ERROR: createDpuKernel(): kernel not yet created");
                return null;
            }
            Expression dim_expr = kernelND.getArgument(2);
            if(dim_expr instanceof IDExpression) {
                dim_expr = Definitions.get((IDExpression)dim_expr);
            }
            
            int dim_kernel = 0;
            if(dim_expr instanceof IntegerLiteral) {
                _kernel.set_num_dim((int)((IntegerLiteral)dim_expr).getValue());
                dim_kernel = (int)((IntegerLiteral)dim_expr).getValue();
            }
            
            CompoundStatement dpu_kernel_call_for_block = new CompoundStatement();
            
            {
            FunctionCall print_call = new FunctionCall(new NameID("printf"), new StringLiteral("Kernel :"+global_kernel_count+" \\n"));
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));
            global_kernel_count++;
            }

            Expression global_expr = kernelND.getArgument(4);
            WorkHorseArray wha = new WorkHorseArray();
            wha.set_ndimension(dim_expr);
            
            if(global_expr instanceof IDExpression) {
                wha.global_work = global_expr;

                FunctionCall print_call = new FunctionCall(new NameID("printf"), new StringLiteral("Global : "));
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));

                for(int dim=0; dim<dim_kernel; dim++) {
                    print_call = new FunctionCall(new NameID("printf"), new StringLiteral("["+dim+"]%ld "), new ArrayAccess(global_expr.clone(), new IntegerLiteral(dim)));
                    dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));
                }
                print_call = new FunctionCall(new NameID("printf"), new StringLiteral("\\n"));
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));

                // wha.set_global_limit(getExpansion((IDExpression)global_expr));

//                 for(Expression e: wha.global_limits) {
//                     System.out.println(e.toString());
//                 }
            }
            
            Expression local_expr = kernelND.getArgument(5);
            if(local_expr instanceof Expression) {
                if(local_expr instanceof Typecast) {
                    local_expr = ((Typecast)local_expr).getExpression();
                }
                
                if(!(local_expr instanceof IDExpression)) {
                    // System.out.println(local_expr.getClass().getName());
                    // for(int arg_ind=0; arg_ind<kernelND.getNumArguments(); arg_ind++) {
                    //     System.out.println(kernelND.getArgument(arg_ind));
                    // }
                    if( (local_expr instanceof IntegerLiteral) && (((IntegerLiteral)local_expr).getValue()==0) ) {
                        _kernel.set_local_dimensions(true);

                        local_expr = new NameID("ss_local_work");

                    }
                }
                FunctionCall print_call = new FunctionCall(new NameID("printf"), new StringLiteral("Local : "));
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));

                for(int dim=0; dim<dim_kernel; dim++) {
                    print_call = new FunctionCall(new NameID("printf"), new StringLiteral("["+dim+"]%ld "), new ArrayAccess(local_expr.clone(), new IntegerLiteral(dim)));
                    dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));
                }
                print_call = new FunctionCall(new NameID("printf"), new StringLiteral("\\n"));
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(print_call));
                

                wha.local_work = local_expr;

                // if(local_expr!=null)
                //     wha.set_local_limit(getExpansion((IDExpression)local_expr));
                // else
                //     wha.set_local_limit(null);

            }
            
            _kernel.set_work_horse_array(wha);

            FunctionCall copy_bd = createDpuCopy(new NameID("_bd"), new IntegerLiteral(0), new NameID("_bd"), null, null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(copy_bd));
            
            FunctionCall func_call = null;
            FunctionCall func_call_offset = null;
            
            CompoundStatement dpu_kernel_call_block = new CompoundStatement(); 
            TranslationUnitKernelInfo t_info = translation_unit_info.get(_kernel.get_kernel_name_string());
            Primitive k_primitive = t_info.primitive;
            
            /*func_call = createDpuCopy(global_expr, new IntegerLiteral(0), global_expr, null, ArgType.NORMAL, DpuCopy.DPU_COPY_TO);
            dpu_kernel_call_block.addStatement(new ExpressionStatement(func_call));
            
            func_call = createDpuCopy(local_expr, new IntegerLiteral(0), local_expr, null, ArgType.NORMAL, DpuCopy.DPU_COPY_TO);
            dpu_kernel_call_block.addStatement(new ExpressionStatement(func_call));*/
            
            //Map<Integer, Pair<Integer, IDExpression>> sorted_map = new TreeMap(_kernel.argument_list);
            
            
            
            Iterator<Map.Entry<Integer, Pair<ArgType, IDExpression>>> arg_iter = _kernel.argument_list.entrySet().iterator();
             
            System.out.println("Kernel Name : " + _kernel.get_kernel_name());

            System.out.println("Iterator size= " + k_primitive.iterator_var.size());
            for(Map.Entry<IDExpression, RangeExpression> iter_var: k_primitive.iterator_var.entrySet()) {
                System.out.println("ITER VARIABLE:" + iter_var.getKey());
                List<Expression> arr_init_val = new ArrayList<Expression>(Arrays.asList(iter_var.getValue().getLB().clone(), iter_var.getValue().getUB().clone()));
                IDExpression arr_iter_name = new NameID("host_"+iter_var.getKey().getName());
                dpu_kernel_call_for_block.addDeclaration(createDecls.createArrayVariableDeclaration(Specifier.INT, arr_iter_name, arr_init_val));

            }

            while(arg_iter.hasNext()) {
                Map.Entry<Integer, Pair<ArgType, IDExpression>> entry = arg_iter.next();
                
                int arg_index = entry.getKey();
                Pair<ArgType, IDExpression> entry_value = entry.getValue();
                func_call = null;
                func_call_offset = null;
                
                IDExpression arg_name = entry_value.getSecond();
                
                switch(entry_value.getFirst()) {
                    case GLOBAL_BUFFER: {
                        // func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), new IntegerLiteral(0), arg_name, null, ArgType.GLOBAL_BUFFER, null, DpuCopy.DPU_COPY_TO);
                        
                        break;
                    }
                        
                    case LOCAL_BUFFER: {
//                         func_call = createDpuCopy(arg_name, new IntegerLiteral(0), arg_name, null, ArgType.NORMAL, DpuCopy.DPU_COPY_TO);
                        break;
                    }
                    case NORMAL: {
                        // IDExpression id_exp = new NameID(t_info.primitive.get_parameter_idexp(arg_index).getName()+"_size");
                        // dpu_kernel_call_for_block.addDeclaration(new VariableDeclaration(Specifier.INT, new VariableDeclarator(id_exp.clone())));
                        // dpu_kernel_call_for_block.addStatement(new ExpressionStatement(new AssignmentExpression(id_exp.clone(), AssignmentOperator.NORMAL, _kernel.get_parameter_size(arg_index).clone())));

                        // FunctionCall loc_func_call = createDpuCopy(id_exp.clone(), new IntegerLiteral(0), id_exp.clone(), null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO);
                        // dpu_kernel_call_for_block.addStatement(makeDpuAssert(loc_func_call));
                        
                        // dpu_kernel_call_for_block.addStatement(createDecls.createCodeAnnotStmt("\n"));
//                         func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), new IntegerLiteral(0), arg_name, )
                    }
                }
                
                if(entry_value.getFirst() == ArgType.GLOBAL_BUFFER) {
                    
//                     System.out.println("entry-key = "+entry.getKey());
                    
                    
                    if(t_info == null) {
                        System.out.println("NULL");
                    }
                    
                    List< Triple<Statement, Integer, AccessType> > offset_start_stmts = _kernel.get_offset_expression(entry.getKey(), arg_name, dpu_kernel_call_for_block, "", true, false); // numder of index types that are gettig used to access this argument
                    List< Triple<Statement, Integer, AccessType> > offset_end_stmts = _kernel.get_offset_expression(entry.getKey(), arg_name, dpu_kernel_call_for_block, "", false, false);
                    
                    

                    if((offset_start_stmts == null) || (offset_start_stmts.size()==0))
                        continue;                    

                    CreateCopyWriteStatements(t_info, _kernel, arg_index, arg_name, null, offset_start_stmts, offset_end_stmts, DpuCopy.DPU_COPY_TO, dpu_kernel_call_for_block);
                }

                else if(entry_value.getFirst() == ArgType.LOCAL_BUFFER) {
                    IDExpression local_mem_dpu_id = new NameID(t_info.primitive.total_parameter_ind_map.get(arg_index).getName()+"_size");
                    Declaration decl = createDecls.createVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), local_mem_dpu_id, _kernel.get_parameter_size(arg_index));
                    dpu_kernel_call_for_block.addDeclaration(decl);
                    func_call = createDpuCopy(local_mem_dpu_id.clone(), new IntegerLiteral(0), local_mem_dpu_id.clone(), null, null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO);
                    dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

                }
            }

            return dpu_kernel_call_for_block;
            
            
        }

        ArrayAccess createAA(IDExpression _id) {
            IDExpression iter_id = new NameID("copy_i");
            return new ArrayAccess(_id.clone(), iter_id.clone());
        }

        ArrayAccess createAAA(int ind) {
            IDExpression se_id = new NameID("ss_SE");
            IDExpression iter_id = new NameID("copy_i");

            return createDecls.createArrayAccess(se_id, iter_id, new IntegerLiteral(ind));
        }

        void CreateCopyWriteStatements(TranslationUnitKernelInfo t_info, KernelRepr _kernel, int arg_index, IDExpression arg_name, Expression cpu_buf, List< Triple<Statement, Integer, AccessType> > offset_start_stmts, List< Triple<Statement, Integer, AccessType> > offset_end_stmts, DpuCopy copy_direction, CompoundStatement dpu_kernel_call_for_block) 
        {
            int num_off_exprs = offset_start_stmts.size();
            List<VariableDeclarator> var_declr = new ArrayList<VariableDeclarator>();

            String pre="";
            if((copy_direction==DpuCopy.DPU_COPY_FROM) || (copy_direction==DpuCopy.DPU_BROADCAST_FROM))
                pre = "p_";

            IndexToAccessInfo k_indInfo = t_info.ind_to_access_info.get(arg_index);
            Primitive k_primitive = t_info.primitive;
                        
            // var_declr.add(new VariableDeclarator(new NameID(pre+arg_name.getName()+"_offset")) );
            // var_declr.add(new VariableDeclarator(new NameID(pre+arg_name.getName()+"_end")) );
            // var_declr.add(new VariableDeclarator(new NameID(pre+arg_name.getName()+"_size")) );

            Literal num_off_exprs_id = new IntegerLiteral(num_off_exprs);

            CompoundStatement copy_comp_stmt = new CompoundStatement();

            IDExpression offset_id = new NameID(pre+arg_name.getName()+"_offset");
            IDExpression offset_cpu_id = new NameID(pre+arg_name.getName()+"_cpu_offset");
            IDExpression end_id = new NameID(pre+arg_name.getName()+"_end");
            IDExpression stride_id = new NameID(pre+arg_name.getName()+"_stride");
            IDExpression final_stride_id = new NameID("f_"+stride_id.getName());
            IDExpression size_id = new NameID(pre+arg_name.getName()+"_size");
            IDExpression map_name = new NameID(pre+arg_name.getName()+"_map");
            IDExpression flag_id = new NameID(pre+arg_name.getName()+"_flag");
            IDExpression length_id = new NameID(pre+arg_name.getName()+"_length");

            IDExpression dpu_arg_name = k_primitive.total_parameter_ind_map.get(arg_index);
            IDExpression dpu_offset_id = new NameID(dpu_arg_name.getName()+"_offset");
            IDExpression dpu_offset_cpu_id = new NameID(dpu_arg_name.getName()+"_cpu_offset");
            IDExpression dpu_stride_id = new NameID(dpu_arg_name.getName()+"_stride");
            IDExpression dpu_map_name = new NameID(dpu_arg_name.getName()+"_map");
            IDExpression dpu_flag_id = new NameID(dpu_arg_name.getName()+"_flag");
            IDExpression dpu_length_id = new NameID(dpu_arg_name.getName()+"_length");



            var_declr.add(new VariableDeclarator(final_stride_id));

            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, offset_id, num_off_exprs, new IntegerLiteral(-1)));
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, end_id, num_off_exprs, new IntegerLiteral(-1)));
            // var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, stride_id, num_off_exprs, new IntegerLiteral(1)));
            VariableDeclarator size_decl = new VariableDeclarator(size_id);
            size_decl.setInitializer(new Initializer(new IntegerLiteral(0)));
            var_declr.add(size_decl);

            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, offset_cpu_id, num_off_exprs, new IntegerLiteral(-1)));
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, length_id, num_off_exprs, new IntegerLiteral(0)));

            // var_declr.add(new VariableDeclarator(, new ArraySpecifier(num_off_exprs_id.clone()) ));
            // var_declr.add(new VariableDeclarator(, new ArraySpecifier(num_off_exprs_id.clone()) ));
            // var_declr.add(new VariableDeclarator(, new ArraySpecifier(num_off_exprs_id.clone()) ));
            // var_declr.add(new VariableDeclarator(, new ArraySpecifier(num_off_exprs_id.clone()) ));

            TranslationUnit k_tran_unit = t_info.tran_unit;
            List<Specifier> dpu_host_specs = new ArrayList<Specifier>(Arrays.asList(new UserSpecifier(new NameID("__host")), new UserSpecifier(new NameID("__dma_aligned")), Specifier.LONG, Specifier.LONG));
            k_tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), dpu_offset_id.clone()));
            k_tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), dpu_offset_cpu_id.clone()));
            // k_tran_unit.addDeclarationFirst(createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), stride_id.clone()));
            k_tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), dpu_flag_id.clone()));
            k_tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), dpu_length_id.clone()));
            k_tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), dpu_stride_id.clone()));
            k_tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createArrayVariableDeclaration(dpu_host_specs, new IntegerLiteral(num_off_exprs), dpu_map_name.clone()));


            Specifier arg_type = k_primitive.total_parameter_spec_map.get(arg_index);

            IDExpression ele_size_id = new NameID(pre+arg_name.getName()+"_ele_size");
            VariableDeclarator ele_size_declr = new VariableDeclarator(ele_size_id);
            ele_size_declr.setInitializer(new Initializer(new SizeofExpression(new ArrayList<Specifier>(Arrays.asList(arg_type)))));
            var_declr.add(ele_size_declr);


            Iterator<Triple<Statement, Integer, AccessType>> offset_start_iterator = offset_start_stmts.iterator(); 
            Iterator<Triple<Statement, Integer, AccessType>> offset_end_iterator = offset_end_stmts.iterator();
            FunctionCall func_call = null;
            
            IDExpression dpu_start_id = new NameID(t_info.primitive.get_parameter_idexp(arg_index).getName()+"_start");
            func_call = createDpuCopy(dpu_start_id, new IntegerLiteral(0), DPU_MRAM_OFFSET.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(makeDpuAssert(func_call));

            Expression cpu_buffer_id = gpu_buffers.get_buffer_map((IDExpression)arg_name);

            // System.out.println("arg_name = " + arg_name);
            if(cpu_buffer_id==null) {
                System.out.println("[CHECK] no matching cpu buffer id for gpu buffer :"+arg_name + " " + gpu_buffers.get_index_global_buffer(arg_name));
                // System.out.println(Long.toBinaryString(_kernel.get_pending_read_status()));
                // long pr = _kernel.get_pending_read_status();
                // for(int _pr=0; _pr<9; _pr++) {
                //     if(((1l<<_pr)&pr) !=0) {
                //         System.out.println("PEN: "+gpu_buffers.get_global_buffer_id(_pr));
                //     }
                // }

                // return;
            }

            // Statement ele_size = createElementSizeStmt(cpu_buffer_id, ele_size_id);
            Statement ele_size = createElementSizeStmt(ele_size_id, k_primitive.total_parameter_spec_map.get(arg_index));
            dpu_kernel_call_for_block.addStatement(ele_size);

            Expression symbol_offset = DPU_MRAM_OFFSET.clone();                                     // new IntegerLiteral(0);
            AssignmentOperator assign_op = AssignmentOperator.NORMAL;

            Long prev_cond = 0l;
            IfStatement prev_if = null;

                        
            Expression[] off_stride = new Expression[num_off_exprs]; 
            Expression[] off_flag = new Expression[num_off_exprs];

            IDExpression temp_buff_id = new NameID("ss_temp"); 
            int W_till = -1;
            int WR_till = -1;


            for(int i=0; i<num_off_exprs; i++) {
                off_stride[i] = new IntegerLiteral(1);

                // start
                Triple<Statement, Integer, AccessType> pair_start = offset_start_iterator.next();
                ExpressionStatement offset_start_expr = (ExpressionStatement)pair_start.getFirst();
                AssignmentExpression start_expr = (AssignmentExpression)offset_start_expr.getExpression();


//                         dpu_kernel_call_for_block.addStatement(_kernel.get_offset_expression_end(arg_name, i, end_expr.getRHS(), ""));       // correct
                //end
                Triple<Statement, Integer, AccessType> pair_end = offset_end_iterator.next();
                ExpressionStatement offset_end_expr = (ExpressionStatement)pair_end.getFirst();
                AssignmentExpression end_expr = (AssignmentExpression)offset_end_expr.getExpression();
                // if(start_expr.getRHS().compareTo(end_expr.getRHS())==0) {
                //     System.out.println("CompareTo (true): " + start_expr);
                //     end_expr.getRHS().swapWith( new BinaryExpression(end_expr.getRHS().clone(), BinaryOperator.ADD, new IntegerLiteral(1)) );
                // }

                Literal flag_val = null;
                switch(pair_start.getThird()) {
                    case WRITE:
                        flag_val = new IntegerLiteral(2);
                        break;
                    case READ:
                        flag_val = new IntegerLiteral(1);
                        break;
                    case WRITE_READ:
                        flag_val = new IntegerLiteral(3);
                        break;
                    default:
                        flag_val = new IntegerLiteral(0);

                }

                off_flag[pair_start.getSecond()] = flag_val;
                off_stride[pair_start.getSecond()] = new IntegerLiteral(1);
                // Expression flag_expr = new AssignmentExpression(flag_id, AssignmentOperator.NORMAL, flag_val);
                // ExpressionStatement offset_flag_expr = new ExpressionStatement(flag_expr);

                ArrayAccess stride_acc = new ArrayAccess(new NameID(pre+arg_name.getName()+"_stride"), new IntegerLiteral(pair_start.getSecond()));
                // Expression stride_expr = new AssignmentExpression(stride_acc, AssignmentOperator.NORMAL, new IntegerLiteral(1));
                // Statement offset_stride_expr = new ExpressionStatement(stride_expr);

                Long cur_cond = k_indInfo.conditional_map.get(pair_start.getSecond());
                // Strides cur_stride = t_info.ind_to_access_info.get(arg_index).stride_info.get(pair_start.getSecond());
                
                WorkFuncRangeMap wf_range_map = new WorkFuncRangeMap();
                wf_range_map.addMapKeys(t_info.conditionals_kernel.convert_ones);

                // if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                    Triple<IfStatement, CompoundStatement, Boolean> insert_real_copy_pair = t_info.conditionals_kernel.constructNestedIfStmt(cur_cond, prev_cond, dpu_kernel_call_for_block, prev_if, wf_range_map, _kernel, k_primitive);

                    System.out.println("WorkFuncRangeMap ("+arg_name+"): ");
                    wf_range_map.print();

                    CompoundStatement insert_real_copy = insert_real_copy_pair.getSecond();

                    if(insert_real_copy_pair.getThird()==true) {
                        List<Declaration> wf_decls = wf_range_map.createWorkFuncDecls(_kernel, k_primitive);
                        for(Declaration new_v_decl: wf_decls) {
                            insert_real_copy.addDeclaration(new_v_decl.clone());
                        }
                    }

                    prev_if = insert_real_copy_pair.getFirst();

                    insert_real_copy.addStatement(offset_start_expr);
                    insert_real_copy.addStatement(offset_end_expr);
                    // insert_real_copy.addStatement(offset_flag_expr);
                    // insert_real_copy.addStatement(offset_stride_expr);
                // }
                  

                


                Expression ana_expr = AnalyseArrayExpression.removeArrayAccess(start_expr.getRHS().clone());
                Strides cur_stride = new Strides();
                cur_stride.findExprStrides(ana_expr, t_info.primitive, _kernel);

                
                Expression size_expr = KernelRepr.get_offset_expression_size(arg_name, ele_size_id, pre, true);    // end-start
                Expression size_expr_byte = new BinaryExpression(size_expr.clone(), BinaryOperator.MULTIPLY, ele_size_id.clone());          // (end-start)*size

                if(i==1) {
                    // symbol_offset = new BinaryExpression(symbol_offset, BinaryOperator.ADD, size_id);
                    symbol_offset = new BinaryExpression(DPU_MRAM_OFFSET.clone(), BinaryOperator.ADD, new BinaryExpression(size_id.clone(), BinaryOperator.MULTIPLY, ele_size_id.clone()));
                    assign_op = AssignmentOperator.ADD;         // dpu_mram_offset+ind_size*ele_ind_size
                }
                
                if((i!=0) && (copy_direction==DpuCopy.DPU_COPY_TO)) {
                    IDExpression size_offset = new NameID(t_info.primitive.get_parameter_idexp(arg_index).getName()+"_offset"+Integer.toString(pair_start.getSecond()));
                    FunctionCall offset_call = createDpuCopy(size_offset, new IntegerLiteral(0), size_id, null, null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO);
                    // insert_real_copy.addStatement(makeDpuAssert(offset_call));
                    // copy_comp_stmt.addStatement(makeDpuAssert(offset_call));
                }



                // if((pair_start.getThird()==AccessType.READ) || (pair_start.getThird()==AccessType.WRITE_READ)) {
                //REAL copy

                if(pair_start.getThird()==AccessType.WRITE)
                // {
                    W_till = i;
                //     WR_till = i;
                // }
                // else if(pair_start.getThird()==AccessType.WRITE_READ) {
                //     WR_till = i;
                // }

                if(copy_direction==DpuCopy.DPU_COPY_TO) {
                    // Expression new_cpu_buffer_id = cpu_buffer_id;
                    func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), symbol_offset, arg_name, null, size_expr_byte, ArgType.GLOBAL_BUFFER, null, copy_direction);

                    if(cur_stride.pattern == 1) {
                        
                        // IDExpression new_stride_name = new NameID(arg_name+"_stride_"+Integer.toString(pair_start.getSecond()));
                        Expression new_stride_name = new ArrayAccess(new NameID(arg_name+"_stride"), new IntegerLiteral(pair_start.getSecond()));

                        // Pair<CompoundStatement, Expression> col_stmt = cur_stride.generateColAccCopy(arg_type, pair_start.getSecond(), cpu_buffer_id, arg_name, temp_buff_id.clone(), copy_direction, func_call, wf_range_map);
                        Pair<Expression, Expression> expr_nstride = cur_stride.createStridedExpr(true);
                        if(expr_nstride==null)
                            System.out.println("[ERROR] new stride is null");
                        
                        // insert_real_copy.addStatement(col_stmt.getFirst());
                        // copy_comp_stmt.addStatement(col_stmt.getFirst());
                        Expression col_stride = cur_stride.stride[cur_stride.pattern_ind[1]];
                        off_stride[pair_start.getSecond()] = col_stride;

                        Expression original_expr = k_indInfo.original_access_map.get(pair_start.getSecond());
                        
                        System.out.println("FIND: "+original_expr);
                        List<Expression> to_modify_exprs = IRTools.findExpressions(k_tran_unit, original_expr);
                        for(Expression m_expr : to_modify_exprs) {
                            ArrayAccess expr_arr_acc = IRTools.getAncestorOfType(m_expr, ArrayAccess.class);
                            if((expr_arr_acc==null)||(!expr_arr_acc.getArrayName().equals(arg_name))) {
                                continue;
                            }
                            if(expr_arr_acc!=null)
                                System.out.println("Array Name: "+expr_arr_acc.getArrayName()+"[] "+arg_name);
                            System.out.println("Before Mod: "+ m_expr);
                            IRTools.replaceAll(m_expr, col_stride, new_stride_name.clone());
                            // Declaration stride_decl = createDecls.createVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, new UserSpecifier(new NameID("__host")), new UserSpecifier(new NameID("__dma_aligned")))), new_stride_name.clone(), new IntegerLiteral(1));
                            // k_tran_unit.addDeclarationFirst(stride_decl);
                            // stride_decl = createDecls.createVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG)), new_stride_name.clone(), new IntegerLiteral(1));
                            
                            // insert_real_copy.addDeclaration(stride_decl.clone());
                            
                            System.out.println("After Mod: "+ m_expr);
                        }

                        // new_cpu_buffer_id = temp_buff_id;

                    }

                    else if(func_call!=null) {
                        // dpu_kernel_call_for_block.addStatement(makeDpuAssert(func_call));
                        // insert_real_copy.addStatement(makeDpuAssert(func_call));
                        // copy_comp_stmt.addStatement(makeDpuAssert(func_call));
                    }

                }
                else {
                    func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), symbol_offset, cpu_buf, null, size_expr_byte, ArgType.GLOBAL_BUFFER, arg_name, copy_direction);
                    if(cur_stride.pattern == 1) {
                        // copy_comp_stmt.addStatement(cur_stride.generateColAccCopy(arg_type, pair_start.getSecond(), cpu_buf, arg_name, new NameID("ss_temp"), copy_direction, func_call, wf_range_map).getFirst());
                        Expression col_stride = cur_stride.stride[cur_stride.pattern_ind[1]];
                        off_stride[i] = col_stride;
                    }
                    else if(func_call!=null) {
                        // dpu_kernel_call_for_block.addStatement(makeDpuAssert(func_call));
                        // copy_comp_stmt.addStatement(makeDpuAssert(func_call));
                    }
                }

                
                // }
                
                // copy_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(size_id.clone(), assign_op, size_expr.clone())));                        
                
                // if((i==num_off_exprs-1) && (copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO)) {
                //     Expression inc_dpu_offset = new BinaryExpression(size_id.clone(), BinaryOperator.MULTIPLY, ele_size_id.clone());
                //     // insert_real_copy.addStatement(new ExpressionStatement(new AssignmentExpression(DPU_MRAM_OFFSET.clone(), AssignmentOperator.ADD, inc_dpu_offset) ));
                //     // insert_real_copy.addStatement(createDecls.createCodeAnnotStmt("\n"));
                //     copy_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(DPU_MRAM_OFFSET.clone(), AssignmentOperator.ADD, inc_dpu_offset) ));
                //     copy_comp_stmt.addStatement(createDecls.createCodeAnnotStmt("\n"));
                // }

                wf_range_map.modify_expr(insert_real_copy);

                prev_cond = cur_cond;
                
            }

            Expression final_stride  = checkUniquenessStrides(off_stride);
            _kernel.setIndexStridePatterns(arg_index, offset_start_stmts.size(), final_stride);

            if(final_stride==null) {
                System.out.println("[ERROR] combination of multiple column strides are present");
                final_stride = new IntegerLiteral(1);
            }
            IDExpression iter_id = new NameID("copy_i");

            Statement ref_stmt = null;
            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                Pair<Statement, Expression> col_stmt = Strides.generateGeneralisedColAccCopy(arg_type, cpu_buffer_id, arg_name, temp_buff_id.clone(), copy_direction, func_call, final_stride, new IntegerLiteral(W_till));
                // copy_comp_stmt.addStatement(col_stmt.getFirst());
                ref_stmt = col_stmt.getFirst();
                // copy_comp_stmt = col_stmt.getFirst();
            }
            else {
                Pair<Statement, Expression> col_stmt = Strides.generateGeneralisedColAccCopy(arg_type, cpu_buf, arg_name, new NameID("ss_temp"), copy_direction, func_call, final_stride, null);
                // copy_comp_stmt.addStatement(col_stmt.getFirst());
                ref_stmt = col_stmt.getFirst();
                // copy_comp_stmt = col_stmt.getFirst();
            }

            Statement copy_stmt = null;

            AssignmentExpression save_start = new AssignmentExpression(new ArrayAccess(offset_cpu_id.clone(), iter_id.clone()), AssignmentOperator.NORMAL, new ArrayAccess(offset_id.clone(), iter_id.clone()));

            Expression off_check = new BinaryExpression(new ArrayAccess(offset_id.clone(), iter_id.clone()), BinaryOperator.COMPARE_LT, new IntegerLiteral(0));
            Statement if_stmt_off_check = new IfStatement(off_check, new ContinueStatement());

            IDExpression se_id = new NameID("ss_SE");
            Expression modify_length_expr = new BinaryExpression(createDecls.createArrayAccess(se_id.clone(), iter_id.clone(), new IntegerLiteral(3)), BinaryOperator.SUBTRACT, createDecls.createArrayAccess(se_id.clone(), iter_id.clone(), new IntegerLiteral(1)));
            modify_length_expr = OneArith.addOne(modify_length_expr);

            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                modify_length_expr = new BinaryExpression(modify_length_expr, BinaryOperator.MULTIPLY, new ArrayAccess(stride_id.clone(), iter_id.clone()));
            }
            // else {
            //     modify_length_expr = _length;
            // }
            IDExpression temp_length_id = new NameID(temp_buff_id.getName()+"_length");
            Expression modify_length = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, modify_length_expr);
            Expression else_modify_length = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, new NameID("LLONG_MAX"));

            Expression length_expr = new ArrayAccess(length_id.clone(), iter_id.clone());
            Statement offset_length_expr = new ExpressionStatement(new AssignmentExpression(length_expr, AssignmentOperator.NORMAL, KernelRepr.get_offset_expression_size(arg_name, ele_size_id, pre, false)));


            Expression stride_assign = new AssignmentExpression(final_stride_id.clone(), AssignmentOperator.NORMAL, final_stride.clone());
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(stride_assign));

            // if(num_off_exprs>1)  {

                FunctionCall print_dpu_data_info = new FunctionCall(new NameID("printf"), new StringLiteral(pre+arg_name.getName()+": %d : %lld %lld L%lld S%lld F%lld M%lld\\n"), iter_id.clone(), createAA(offset_id), createAA(offset_cpu_id), createAA(length_id), createAA(stride_id), createAA(flag_id), createAA(map_name));
                FunctionCall print_index_map_info = new FunctionCall(new NameID("printf"), new StringLiteral(pre+arg_name.getName()+"_index[%d] : %lld %lld, stride=%lld, SE: %d %d %d %d\\n"), iter_id.clone(), createAA(offset_id), createAA(end_id), final_stride_id.clone(), createAAA(0), createAAA(1), createAAA(2), createAAA(3));

                VariableDeclarator map_declr = new VariableDeclarator(map_name.clone(), ArraySpecifier.UNBOUNDED);
                List<Expression> map_list = new ArrayList<Expression>();
                for(int j=0; j<num_off_exprs; j++) {
                    map_list.add(new IntegerLiteral(j));
                }
                Initializer map_init = new Initializer(map_list);
                map_declr.setInitializer(map_init);

                dpu_kernel_call_for_block.addDeclaration(new VariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), map_declr));
                
                if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                    FunctionCall gen_rect = new FunctionCall(new NameID("generateRectangle"), offset_id.clone(), end_id.clone(), stride_id.clone(), map_name.clone(), flag_id.clone(), final_stride_id.clone(), new IntegerLiteral(num_off_exprs));
                    Statement gen_rect_stmt = createCodeAnnotStmt("ss_SE = (int(*)[4])"+gen_rect+";\n");
                    dpu_kernel_call_for_block.addStatement(gen_rect_stmt);
                }
                

                Expression if_cond = new BinaryExpression(new ArrayAccess(map_name.clone(), iter_id.clone()), BinaryOperator.COMPARE_EQ, iter_id.clone());
                Expression if_cond_2 = new BinaryExpression(se_id.clone(), BinaryOperator.COMPARE_NE, new NameID("NULL"));
                // if_cond = new BinaryExpression(if_cond, BinaryOperator.LOGICAL_AND, if_cond_2);

                CompoundStatement if_path = new CompoundStatement();

                IfStatement modify_length_if = new IfStatement(if_cond_2, new ExpressionStatement(modify_length), new ExpressionStatement(else_modify_length));

                if_path.addStatement(modify_length_if);

                // System.out.println("Ref: "+ref_stmt);
                if_path.addStatement(ref_stmt);
                if_path.addStatement(offset_length_expr);
                if_path.addStatement(createCopyOffsetStride(pre, arg_name, true));

                // Expression size_expr = KernelRepr.get_offset_expression_size(arg_name, ele_size_id, pre, true);
                Expression size_expr = KernelRepr.make8ByteAligned(createAA(length_id), ele_size_id.clone());
                if_path.addStatement(new ExpressionStatement(new AssignmentExpression(size_id.clone(), assign_op, size_expr.clone()))); 
                

                copy_stmt = new IfStatement(if_cond, if_path, createCopyOffsetStride(pre, arg_name, false));
                copy_comp_stmt.addStatement(copy_stmt);
                copy_comp_stmt.addStatementBefore(copy_stmt, if_stmt_off_check);
                // copy_comp_stmt.addStatement(new ExpressionStatement(print_dpu_data_info));
                // copy_comp_stmt.addStatementBefore(copy_stmt, new ExpressionStatement(print_index_map_info));

                // copy_comp_stmt.addStatementBefore(copy_stmt, new ExpressionStatement(modify_length));
            // }

            // else {
            //     // copy_stmt = copy_comp_stmt;
            //     copy_comp_stmt.addStatement(ref_stmt);
            //     copy_comp_stmt.addStatement(offset_length_expr);
            //     copy_comp_stmt.addStatement(createCopyOffsetStride(pre, arg_name, true));
            //     copy_comp_stmt.addStatementBefore(ref_stmt, if_stmt_off_check);
            //     copy_comp_stmt.addStatement(new ExpressionStatement(print_dpu_data_info));
            //     copy_comp_stmt.addStatementBefore(ref_stmt, new ExpressionStatement(print_index_map_info));
            //     copy_comp_stmt.addStatementBefore(ref_stmt, new ExpressionStatement(modify_length));

            //     Expression size_expr = KernelRepr.get_offset_expression_size(arg_name, ele_size_id, pre, true);
            //     copy_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(size_id.clone(), assign_op, size_expr.clone()))); 
                
            // }
            copy_comp_stmt.addStatementBefore(if_stmt_off_check, new ExpressionStatement(save_start));

            dpu_kernel_call_for_block.addDeclaration(new VariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), var_declr));
            dpu_kernel_call_for_block.addDeclaration(createDecls.createArrayVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), stride_id.clone(), new ArrayList<Expression>(Arrays.asList(off_stride)) ));
            dpu_kernel_call_for_block.addDeclaration(createDecls.createArrayVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), flag_id.clone(), new ArrayList<Expression>(Arrays.asList(off_flag)) ));

                
            
            

            
            
            // IDExpression prvs_dpu_mram_offset = new NameID("prvs_"+DPU_MRAM_OFFSET.getName());
            // Expression store_prvs = new AssignmentExpression(prvs_dpu_mram_offset, AssignmentOperator.NORMAL, DPU_MRAM_OFFSET.clone());

            // dpu_kernel_call_for_block.addStatement(new ExpressionStatement(store_prvs));

            Literal limit_val = new IntegerLiteral(num_off_exprs);
            Expression condition = new BinaryExpression(iter_id.clone(), BinaryOperator.COMPARE_LT, limit_val);
            Expression step = new UnaryExpression(UnaryOperator.POST_INCREMENT, iter_id.clone());
            // Statement init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0)));
            Statement init_stmt = new DeclarationStatement(createDecls.createVariableDeclaration(Specifier.INT, iter_id.clone(), new IntegerLiteral(0)));
            ForLoop copy_for_loop = new ForLoop(init_stmt, condition, step, copy_comp_stmt);

            
            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO)
                dpu_kernel_call_for_block.addStatement(copy_for_loop);
            // dpu_kernel_call_for_block.addStatementBefore(copy_for_loop, if_stmt_off_check);

            func_call = createDpuCopy(dpu_offset_id.clone(), new IntegerLiteral(0), offset_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_stride_id.clone(), new IntegerLiteral(0), stride_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_flag_id.clone(), new IntegerLiteral(0), flag_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_length_id.clone(), new IntegerLiteral(0), length_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_offset_cpu_id.clone(), new IntegerLiteral(0), offset_cpu_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            // if(num_off_exprs>1) {
                func_call = createDpuCopy(dpu_map_name.clone(), new IntegerLiteral(0), map_name.clone(), null, null, ArgType.NORMAL, null, copy_direction);
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));
            // }

            if(copy_direction==DpuCopy.DPU_COPY_FROM || copy_direction==DpuCopy.DPU_BROADCAST_FROM)
                dpu_kernel_call_for_block.addStatement(copy_for_loop);

            // if((i==num_off_exprs-1) && (copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO)) {

            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                Expression inc_dpu_offset = new BinaryExpression(size_id.clone(), BinaryOperator.MULTIPLY, ele_size_id.clone());
                // insert_real_copy.addStatement(new ExpressionStatement(new AssignmentExpression(DPU_MRAM_OFFSET.clone(), AssignmentOperator.ADD, inc_dpu_offset) ));
                // insert_real_copy.addStatement(createDecls.createCodeAnnotStmt("\n"));
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(new AssignmentExpression(DPU_MRAM_OFFSET.clone(), AssignmentOperator.ADD, inc_dpu_offset) ));
                dpu_kernel_call_for_block.addStatement(createDecls.createCodeAnnotStmt("\n"));  
            }
                                         

                // }

        }

void CreateCopyReadStatements(TranslationUnitKernelInfo t_info, KernelRepr _kernel, int arg_index, IDExpression arg_name, Expression cpu_buf, DpuCopy copy_direction, CompoundStatement dpu_kernel_call_for_block) 
        {
            int num_off_exprs = _kernel.getNumIndexPatterns(arg_index);

            List<VariableDeclarator> var_declr = new ArrayList<VariableDeclarator>();

            String pre="";
            if((copy_direction==DpuCopy.DPU_COPY_FROM) || (copy_direction==DpuCopy.DPU_BROADCAST_FROM))
                pre = "p_";
                        
            IndexToAccessInfo k_indInfo = t_info.ind_to_access_info.get(arg_index);
            Primitive k_primitive = t_info.primitive;
            // var_declr.add(new VariableDeclarator(new NameID(pre+arg_name.getName()+"_offset")) );
            // var_declr.add(new VariableDeclarator(new NameID(pre+arg_name.getName()+"_end")) );
            // var_declr.add(new VariableDeclarator(new NameID(pre+arg_name.getName()+"_size")) );

            Literal num_off_exprs_id = new IntegerLiteral(num_off_exprs);

            CompoundStatement copy_comp_stmt = new CompoundStatement();

            IDExpression offset_id = new NameID(pre+arg_name.getName()+"_offset");
            IDExpression offset_cpu_id = new NameID(pre+arg_name.getName()+"_cpu_offset");
            IDExpression end_id = new NameID(pre+arg_name.getName()+"_end");
            IDExpression stride_id = new NameID(pre+arg_name.getName()+"_stride");
            IDExpression final_stride_id = new NameID("f_"+stride_id.getName());
            IDExpression size_id = new NameID(pre+arg_name.getName()+"_size");
            IDExpression map_name = new NameID(pre+arg_name.getName()+"_map");
            IDExpression flag_id = new NameID(pre+arg_name.getName()+"_flag");
            IDExpression length_id = new NameID(pre+arg_name.getName()+"_length");

            IDExpression dpu_arg_name = k_primitive.total_parameter_ind_map.get(arg_index);
            IDExpression dpu_offset_id = new NameID(dpu_arg_name.getName()+"_offset");
            IDExpression dpu_offset_cpu_id = new NameID(dpu_arg_name.getName()+"_cpu_offset");
            IDExpression dpu_stride_id = new NameID(dpu_arg_name.getName()+"_stride");
            IDExpression dpu_map_name = new NameID(dpu_arg_name.getName()+"_map");
            IDExpression dpu_flag_id = new NameID(dpu_arg_name.getName()+"_flag");
            IDExpression dpu_length_id = new NameID(dpu_arg_name.getName()+"_length");

            var_declr.add(new VariableDeclarator(final_stride_id));

            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, offset_id, num_off_exprs, null));
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, length_id, num_off_exprs, null));
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, stride_id, num_off_exprs, null));
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, offset_cpu_id, num_off_exprs, null));
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, flag_id, num_off_exprs, null));
            // if(num_off_exprs>1)
            var_declr.add(createDecls.createArrayVariableDeclarator(Specifier.INT, map_name, num_off_exprs, null));

            Specifier arg_type = k_primitive.total_parameter_spec_map.get(arg_index);

            IDExpression ele_size_id = new NameID(pre+arg_name.getName()+"_ele_size");
            // var_declr.add(new VariableDeclarator(ele_size_id));
            VariableDeclarator ele_size_declr = new VariableDeclarator(ele_size_id);
            ele_size_declr.setInitializer(new Initializer(new SizeofExpression(new ArrayList<Specifier>(Arrays.asList(arg_type)))));
            var_declr.add(ele_size_declr);

            FunctionCall func_call = null;
            
            IDExpression dpu_start_id = new NameID(t_info.primitive.get_parameter_idexp(arg_index).getName()+"_start");
            func_call = createDpuCopy(dpu_start_id, new IntegerLiteral(0), DPU_MRAM_OFFSET.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(makeDpuAssert(func_call));

            Expression cpu_buffer_id = gpu_buffers.get_buffer_map((IDExpression)arg_name);

            // Statement ele_size = createElementSizeStmt(cpu_buffer_id, ele_size_id);
            Statement ele_size = createElementSizeStmt(ele_size_id, k_primitive.total_parameter_spec_map.get(arg_index));

            dpu_kernel_call_for_block.addStatement(ele_size);

            Expression size_expr = KernelRepr.get_offset_expression_size(arg_name, ele_size_id, pre, true);    // end-start
            Expression size_expr_byte = new BinaryExpression(size_expr.clone(), BinaryOperator.MULTIPLY, ele_size_id.clone());          // (end-start)*size

            Expression symbol_offset = DPU_MRAM_OFFSET.clone();    

            AssignmentOperator assign_op = AssignmentOperator.NORMAL;

            Long prev_cond = 0l;
            IfStatement prev_if = null;

            
            

            IDExpression temp_buff_id = new NameID("ss_temp"); 
            int W_till = -1;
            int WR_till = -1;

            func_call = createDpuCopy(dpu_offset_id.clone(), new IntegerLiteral(0), offset_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_stride_id.clone(), new IntegerLiteral(0), stride_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_flag_id.clone(), new IntegerLiteral(0), flag_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_length_id.clone(), new IntegerLiteral(0), length_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            func_call = createDpuCopy(dpu_offset_cpu_id.clone(), new IntegerLiteral(0), offset_cpu_id.clone(), null, null, ArgType.NORMAL, null, copy_direction);
            dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));

            // if(num_off_exprs>1) {
                func_call = createDpuCopy(dpu_map_name.clone(), new IntegerLiteral(0), map_name.clone(), null, null, ArgType.NORMAL, null, copy_direction);
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(func_call));
            // }

            Expression final_stride  = _kernel.getFinalStride(arg_index);

            if(final_stride==null) {
                System.out.println("[ERROR] combination of multiple column strides are present");
                final_stride=new IntegerLiteral(1);
            }
            else
                final_stride = _kernel.getFinalStride(arg_index).clone();
            
            IDExpression iter_id = new NameID("copy_i");

            Statement ref_stmt = null;


            symbol_offset = new BinaryExpression(symbol_offset, BinaryOperator.ADD, new BinaryExpression(new ArrayAccess(offset_id.clone(), iter_id.clone()), BinaryOperator.MULTIPLY, ele_size_id.clone()));
            Expression cpu_offset = new ArrayAccess(offset_cpu_id.clone(), iter_id.clone());

            Expression copy_size = KernelRepr.make8ByteAligned(createAA(length_id), ele_size_id.clone());
            copy_size = new BinaryExpression(copy_size, BinaryOperator.MULTIPLY, ele_size_id.clone());

            func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), symbol_offset, cpu_buf, cpu_offset, copy_size, ArgType.GLOBAL_BUFFER, arg_name, copy_direction);
            Pair<Statement, Expression> col_stmt = null;
                col_stmt = Strides.generateGeneralisedColAccCopy(arg_type, cpu_buf, arg_name, new NameID("ss_temp"), copy_direction, func_call, final_stride, null);
                ref_stmt = col_stmt.getFirst();
            // copy_comp_stmt.addStatement(col_stmt.getFirst());
            
            // copy_comp_stmt = col_stmt.getFirst();


            Statement copy_stmt = null;

            Expression off_check = new BinaryExpression(new ArrayAccess(offset_id.clone(), iter_id.clone()), BinaryOperator.COMPARE_LT, new IntegerLiteral(0));
            Statement if_stmt_off_check = new IfStatement(off_check, new ContinueStatement());

            Expression length_expr = new ArrayAccess(length_id.clone(), iter_id.clone());
            Statement offset_length_expr = new ExpressionStatement(new AssignmentExpression(length_expr, AssignmentOperator.NORMAL, KernelRepr.get_offset_expression_size(arg_name, ele_size_id, pre, false)));

            IDExpression temp_length_id = new NameID(temp_buff_id.getName()+"_length");
            Expression modify_length = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, new ArrayAccess(length_id.clone(), iter_id.clone()));
            
            FunctionCall print_dpu_data_info = new FunctionCall(new NameID("printf"), new StringLiteral(pre+arg_name.getName()+": %d : %ld %lld %lld %lld %lld\\n"), iter_id.clone(), createAA(offset_id), createAA(offset_cpu_id), createAA(length_id), createAA(stride_id), createAA(flag_id));


            // if(num_off_exprs>1)  {
                print_dpu_data_info = new FunctionCall(new NameID("printf"), new StringLiteral(pre+arg_name.getName()+": %d : %lld %lld L%lld S%lld F%lld M%lld\\n"), iter_id.clone(), createAA(offset_id), createAA(offset_cpu_id), createAA(length_id), createAA(stride_id), createAA(flag_id), createAA(map_name));

                // VariableDeclarator map_declr = new VariableDeclarator(map_name.clone(), ArraySpecifier.UNBOUNDED);
                // List<Expression> map_list = new ArrayList<Expression>();
                // for(int j=0; j<num_off_exprs; j++) {
                //     map_list.add(new IntegerLiteral(j));
                // }
                // Initializer map_init = new Initializer(map_list);
                // map_declr.setInitializer(map_init);

                // dpu_kernel_call_for_block.addDeclaration(new VariableDeclaration(Specifier.INT, map_declr));
                
                Expression stride_assign = new AssignmentExpression(final_stride_id.clone(), AssignmentOperator.NORMAL, final_stride.clone());
                dpu_kernel_call_for_block.addStatement(new ExpressionStatement(stride_assign));

                if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                    FunctionCall gen_rect = new FunctionCall(new NameID("generateRectangle"), final_stride_id.clone(), new IntegerLiteral(num_off_exprs), offset_id.clone(), end_id.clone(), stride_id.clone(), map_name.clone(), flag_id.clone());
                    Statement gen_rect_stmt = createCodeAnnotStmt("ss_SE = (int(*)[4])"+gen_rect+";\n");
                    dpu_kernel_call_for_block.addStatement(gen_rect_stmt);
                }
                

                Expression if_cond = new BinaryExpression(new ArrayAccess(map_name.clone(), iter_id.clone()), BinaryOperator.COMPARE_EQ, iter_id.clone());
                
                CompoundStatement if_path = new CompoundStatement();
                if_path.addStatement(ref_stmt);

                copy_stmt = new IfStatement(if_cond, if_path);
                copy_comp_stmt.addStatement(copy_stmt);
                copy_comp_stmt.addStatementBefore(copy_stmt, if_stmt_off_check);
                // copy_comp_stmt.addStatementBefore(copy_stmt, new ExpressionStatement(print_dpu_data_info));
                copy_comp_stmt.addStatementBefore(copy_stmt, new ExpressionStatement(modify_length));

            // }

            // else {
            //     // copy_stmt = copy_comp_stmt;
            //     copy_comp_stmt.addStatement(ref_stmt);
            //     copy_comp_stmt.addStatementBefore(ref_stmt, if_stmt_off_check);
            //     copy_comp_stmt.addStatementBefore(ref_stmt, new ExpressionStatement(print_dpu_data_info));
            //     copy_comp_stmt.addStatementBefore(ref_stmt, new ExpressionStatement(modify_length));

            // }


            dpu_kernel_call_for_block.addDeclaration(new VariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), var_declr));

            Literal limit_val = new IntegerLiteral(num_off_exprs);
            Expression condition = new BinaryExpression(iter_id.clone(), BinaryOperator.COMPARE_LT, limit_val);
            Expression step = new UnaryExpression(UnaryOperator.POST_INCREMENT, iter_id.clone());
            // Statement init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0)));
            Statement init_stmt = new DeclarationStatement(createDecls.createVariableDeclaration(Specifier.INT, iter_id.clone(), new IntegerLiteral(0)));

            ForLoop copy_for_loop = new ForLoop(init_stmt, condition, step, copy_comp_stmt);
            
            dpu_kernel_call_for_block.addStatement(copy_for_loop);


        }



        Statement createCopyOffsetStride(String pre, IDExpression arg_name, boolean is_if) {
            IDExpression offset_id = new NameID(pre+arg_name.getName()+"_offset");
            IDExpression end_id = new NameID(pre+arg_name.getName()+"_end");
            IDExpression stride_id = new NameID(pre+arg_name.getName()+"_stride");
            IDExpression final_stride_id = new NameID("f_"+stride_id.getName());
            IDExpression size_id = new NameID(pre+arg_name.getName()+"_size");
            IDExpression map_name = new NameID(pre+arg_name.getName()+"_map");
            IDExpression off_ind_name = new NameID(pre+arg_name.getName()+"_off_ind");
            IDExpression prvs_dpu_mram_offset = new NameID("prvs_"+DPU_MRAM_OFFSET.getName());

            IDExpression copy_id = new NameID("copy_i");
            IDExpression se_id = new NameID("ss_SE");

            Expression offset_ind = new ArrayAccess(off_ind_name.clone(), copy_id.clone());
            
            ArrayAccess offset_acc = new ArrayAccess(offset_id.clone(), copy_id.clone());
            Expression offset_acc_val = null;
            ArrayAccess stride_acc = new ArrayAccess(stride_id.clone(), copy_id.clone());
            Expression stride_acc_val = null;

            Expression stride_assign = null;

            CompoundStatement c_st = new CompoundStatement();

            if(is_if) {
                offset_acc_val = size_id.clone();
            }

            else {

                Expression offset_acc_map = new ArrayAccess(map_name.clone(), copy_id.clone());
                offset_acc_map = new ArrayAccess(offset_id.clone(), offset_acc_map);

                // stride_acc = new ArrayAccess(stride_id.clone(), copy_id.clone());
                stride_acc_val = new ArrayAccess(map_name.clone(), copy_id.clone());
                stride_acc_val = new ArrayAccess(stride_id.clone(), stride_acc_val); 

                offset_acc_val = new BinaryExpression(stride_acc.clone(), BinaryOperator.MULTIPLY, createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(1)));
                offset_acc_val = new BinaryExpression(offset_acc_val, BinaryOperator.ADD, createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(0)));
                offset_acc_val = new BinaryExpression(offset_acc_map.clone(), BinaryOperator.ADD, offset_acc_val);

                stride_assign = new AssignmentExpression(stride_acc, AssignmentOperator.NORMAL, stride_acc_val);
                c_st.addStatement(new ExpressionStatement(stride_assign));

                // offset_acc = new ArrayAccess(offset_id.clone(), copy_id.clone());
            }

            
            Expression offset_assign = new AssignmentExpression(offset_acc, AssignmentOperator.NORMAL, offset_acc_val);
            
            if(is_if)
                return new ExpressionStatement(offset_assign);
            c_st.addStatement(new ExpressionStatement(offset_assign));
            return c_st;

        }


        static Expression createSizeofRectange() {
            IDExpression copy_id = new NameID("copy_i");
            IDExpression se_id = new NameID("ss_SE");

            Expression x_dist = createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(2));
            x_dist = new BinaryExpression(x_dist, BinaryOperator.SUBTRACT, createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(0)));
            x_dist = OneArith.addOne(x_dist);

            Expression y_dist = createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(3));
            y_dist = new BinaryExpression(y_dist, BinaryOperator.SUBTRACT, createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(1)));
            y_dist = OneArith.addOne(y_dist);

            Expression size_expr = new BinaryExpression(x_dist , BinaryOperator.MULTIPLY, y_dist);
            return size_expr;

        }

        Expression checkUniquenessStrides(Expression[] off_stride) {
            Expression return_stride = new IntegerLiteral(1);
            Expression one = new IntegerLiteral(1);
            for(int i=0; i<off_stride.length; i++) {
                if(!off_stride[i].equals(one) && !off_stride[i].equals(return_stride)) {
                    if(return_stride.equals(one))
                        return_stride = off_stride[i];
                    else {
                        System.out.println("!ERROR "+return_stride +" " + off_stride[i]);
                        return null;
                    }
                }
            }

            return return_stride;

        }
        
        CompoundStatement createDpuKernelInit(KernelRepr _kernel) {
        
            if(_kernel==null) {
                System.out.println("ERROR: createDpuKernelInit(): kernel not yet created");
                return null;
            }
            CompoundStatement dpu_kernel_call_block = new CompoundStatement(); 
            // DpuCopy init_transfer_type = DpuCopy.DPU_BROADCAST_TO;
            DpuCopy init_transfer_type = DpuCopy.DPU_COPY_TO;

            FunctionCall func_call = null;
            
            func_call = createDpuCopy(new NameID("_gs"), new IntegerLiteral(0), _kernel.wha.global_work, null, null, ArgType.NORMAL, null, init_transfer_type);
            dpu_kernel_call_block.addStatement(makeDpuAssert(func_call));
            
            
            func_call = createDpuCopy(new NameID("_bs"), new IntegerLiteral(0), _kernel.wha.local_work, null, null, ArgType.NORMAL, null, init_transfer_type);
            dpu_kernel_call_block.addStatement(makeDpuAssert(func_call));
            
            Iterator<Map.Entry<Integer, Pair<ArgType, IDExpression>>> arg_iter = _kernel.argument_list.entrySet().iterator();
            // CompoundStatement dpu_kernel_call_for_block = new CompoundStatement(); 
            
            while(arg_iter.hasNext()) {
                Map.Entry<Integer, Pair<ArgType, IDExpression>> entry = arg_iter.next();
                Pair<ArgType, IDExpression> entry_value = entry.getValue();
                func_call = null;
                IDExpression arg_name = entry_value.getSecond();
                
                switch(entry_value.getFirst()) {
                    case NORMAL: {
                        func_call = createDpuCopy(arg_name, new IntegerLiteral(0), arg_name, null, null, ArgType.NORMAL, null, init_transfer_type);
                        break;
                    }
                }
                if(func_call!=null)
                    dpu_kernel_call_block.addStatement(makeDpuAssert(func_call));
            }
            
            String code = String.format("int N_TASKLETS = %d;\n", 0);
                    
            IDExpression id_exp = null;
            
            
            
            func_call = createDpuCopy(new NameID("dpu_gs"), new BinaryExpression(new IntegerLiteral(_kernel.partition_dim_grid), BinaryOperator.MULTIPLY, createDecls.createSizeofExpression(new NameID("long")) ), new NameID("N_MULTI_WGS"), null, null, ArgType.NORMAL, null, init_transfer_type);
            dpu_kernel_call_block.addStatement(makeDpuAssert(func_call));
            
            
            return dpu_kernel_call_block;
            
            
        }

        CompoundStatement createDpuParamTrace() {
            CompoundStatement trace_params = new CompoundStatement();
            IDExpression func_id = new NameID("note_down");
            
            trace_params.addStatement(new ExpressionStatement(new FunctionCall(func_id.clone(), new NameID("N_TASKLETS")) ));
            trace_params.addStatement(new ExpressionStatement(new FunctionCall(func_id.clone(), new NameID("N_MULTI_WGS")) ));
            trace_params.addStatement(new ExpressionStatement(new FunctionCall(func_id.clone(), new NameID("PARTITION_DIM_GRID")) ));
            trace_params.addStatement(new ExpressionStatement(new FunctionCall(func_id.clone(), new NameID("PARTITION_DIM_WG")) ));

            return trace_params;
        }

        void createDpuParamCopy(KernelRepr _kernel, CompoundStatement block_var, CompoundStatement block_copy)
        {
            createDeclCopy(block_var, block_copy, new NameID("N_TASKLETS"), _kernel.n_tasklets);
            // createDeclCopy(block_var, block_copy, new NameID("MULTI_WGS"), _kernel.multi_wgs);
            
            createDeclCopy(block_var, block_copy, new NameID("PARTITION_DIM_GRID"), _kernel.partition_dim_grid);
            createDeclCopy(block_var, block_copy, new NameID("PARTITION_DIM_WG"), _kernel.partition_dim_wg);


        }

        void createDpuParamDeclarations(KernelRepr _kernel, CompoundStatement ndrange_block, int[] order) {
//             createDeclCopy(ndrange_block, new NameID("N_MULTI_WGS"), _kernel.n_multi_wgs);
            ndrange_block.addDeclaration(createDecls.createVariableDeclaration(Specifier.LONG, new NameID("N_MULTI_WGS"), new IntegerLiteral(_kernel.n_multi_wgs)));
            ndrange_block.addDeclaration(createDecls.createVariableDeclaration(Specifier.LONG, new NameID("Work_Group_Id"), new IntegerLiteral(0)));
            ndrange_block.addDeclaration(createDecls.createVariableDeclaration(Specifier.LONG, new NameID("p_Work_Group_Id"), new IntegerLiteral(0)));

            if(_kernel.local_manual) {
                IDExpression local_expr = new NameID("ss_local_work");
                        
                List<Expression> local_dimensions = new ArrayList<Expression>();
                for(int dim=0; dim<_kernel.num_dim; dim++) {
                    local_dimensions.add(new IntegerLiteral(_kernel.manual_local_dimensions[dim]));
                }

                Declaration decl = createDecls.createArrayVariableDeclaration(Specifier.LONG, (IDExpression)local_expr.clone(), local_dimensions);
                System.out.println("LOCAL_DECL: " + decl);
                ndrange_block.addDeclaration(decl);    
            }
            
            VariableDeclarator decl_bd = null;
            VariableDeclaration dec_b = null;
            
            for(int dim=0; dim<3; dim++) {
                decl_bd = new VariableDeclarator(getIteratorVariable(PolyStart.B_START, dim, ""));
                decl_bd.setInitializer(new Initializer(new IntegerLiteral(0)));
                if(dim==0)
                    dec_b = new VariableDeclaration(Specifier.LONG, decl_bd);
                else
                    dec_b.addDeclarator(decl_bd);
                
                decl_bd = new VariableDeclarator(getIteratorVariable(PolyStart.B_START, dim, "p_",""));
                decl_bd.setInitializer(new Initializer(new IntegerLiteral(0)));
                dec_b.addDeclarator(decl_bd);
            }

            ndrange_block.addDeclaration(dec_b);
            
            List<Statement> outer_stmts = createBlockLimits(_kernel.wha, order, true);
            for(Statement st: outer_stmts) {
                if(st instanceof DeclarationStatement) {
                    ndrange_block.addDeclaration(((DeclarationStatement)st).getDeclaration().clone());
                }
                else
                    ndrange_block.addStatement(st);
            }
            
            TranslationUnitKernelInfo t_info = translation_unit_info.get(_kernel.get_kernel_name_string());
            if(!t_info.is_initialized()) {
                String code = "";
                int div = _kernel.n_tasklets/_kernel.n_multi_wgs;
                int rem = _kernel.n_tasklets%_kernel.n_multi_wgs;
                int _div = div+1;
//                 int rem = _kernel.wha.global_limits_int[_kernel.partition_dim_grid]%_kernel.n_tasklets;
                String code_array = "barrier_t* barrier[] = {";
                for(int ind=0; ind<_kernel.n_multi_wgs; ind++) {
                    if(ind < rem) {
                        code += "BARRIER_INIT(barrier_"+ind+", "+ _div +");\n";
                    }
                    else {
                        code += "BARRIER_INIT(barrier_"+ind+", "+ div +");\n";
                    }
                    if(ind>0)
                        code_array += ", ";
                    code_array += "&barrier_"+ind;
                }
                code_array += "};\n";
                String code_define = "#define N_MULTI_WGS "+_kernel.n_multi_wgs +"\n";
                
                t_info.tran_unit.addDeclarationAfter(t_info.getRefDeclaration(), createDecls.createCodeAnnotDecl(code+code_array+code_define));
//                 t_info.done_intialization(true);
            }
            
        }
        
        
        
        void createDeclCopy(CompoundStatement c_stmt_var, CompoundStatement c_stmt_copy, IDExpression id_exp, int val) {
            c_stmt_var.addDeclaration(createDecls.createVariableDeclaration(Specifier.LONG, id_exp, new IntegerLiteral(val)));
            c_stmt_copy.addStatement(makeDpuAssert(createDpuCopy(id_exp.clone(), new IntegerLiteral(0), id_exp.clone(), null, null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO)));
        }
        
        static FunctionCall createDpuCopy(Expression dpu_sym_id, Expression dpu_sym_offset, Expression cpu_sym_id, Expression cpu_offset, Expression size, ArgType is_pointer, Expression gpu_buf_in_cpu, DpuCopy from) {
            // System.out.println("Sym Ids "+dpu_sym_id+" "+cpu_sym_id);
            FunctionCall fun_call = null;
            if(from==DpuCopy.DPU_COPY_FROM)
                fun_call = new FunctionCall(new NameID("dpu_copy_from"));
            else if(from==DpuCopy.DPU_COPY_TO)
                fun_call = new FunctionCall(new NameID("dpu_copy_to"));
            else if(from==DpuCopy.DPU_BROADCAST_FROM)
                fun_call = new FunctionCall(new NameID("dpu_broadcast_from"));
            else
                fun_call = new FunctionCall(new NameID("dpu_broadcast_to"));
                
            fun_call.addArgument(new NameID(dpu_name));
            if(is_pointer==ArgType.GLOBAL_BUFFER){
                fun_call.addArgument(new NameID("DPU_MRAM_HEAP_POINTER_NAME"));
                // fun_call.addArgument(new StringLiteral("GLOBAL_BUFFER"));
                // fun_call.addArgument(new NameID("dpu_mram_offset"));
            }
            else {
                fun_call.addArgument(new StringLiteral(dpu_sym_id.toString()));
            }
            fun_call.addArgument(dpu_sym_offset.clone());
            if(is_pointer==ArgType.GLOBAL_BUFFER) {
                BinaryExpression bin_expr = null;
//                 if(offset_ind == 0)
        
                    Expression glob_cpu_sym_id = cpu_sym_id;
                    
                    if((from == DpuCopy.DPU_COPY_TO) || (from == DpuCopy.DPU_BROADCAST_TO)) {
                        if(cpu_sym_id instanceof IDExpression) {
                            Expression temp_sym_id= gpu_buffers.get_buffer_map((IDExpression)cpu_sym_id);
                            if(temp_sym_id!=null)
                                glob_cpu_sym_id = temp_sym_id;
                            // else if(from == DpuCopy.DPU_COPY_TO)
                            else
                                return null;
                        }
                    }
                    
                    
                    
                    Expression offset_expr = cpu_offset;
                    Expression copy_id = new NameID("copy_i");
                    if(offset_expr==null) {
                        if((from == DpuCopy.DPU_COPY_FROM) || (from == DpuCopy.DPU_BROADCAST_FROM))
                            offset_expr = new ArrayAccess(new NameID("p_"+gpu_buf_in_cpu.toString()+"_offset"), copy_id.clone());
                        else
                            offset_expr = new ArrayAccess(new NameID(cpu_sym_id.toString()+"_offset"), copy_id.clone());
                    }
                    bin_expr = new BinaryExpression(glob_cpu_sym_id.clone(), BinaryOperator.ADD, offset_expr);
//                 else
//                     bin_expr = new BinaryExpression(cpu_sym_id.clone(), BinaryOperator.ADD, new NameID(cpu_sym_id.toString()+"_offset_"+Integer.toString(offset_ind)));
                    
                fun_call.addArgument(bin_expr);
                if(size==null)
                    fun_call.addArgument(new NameID(cpu_sym_id.toString()+"_partition_size"));
                else
                    fun_call.addArgument(size);
            }
            else if(is_pointer==ArgType.NORMAL || is_pointer==ArgType.LOCAL_BUFFER){
                if(size==null) {
                    if(is_pointer == ArgType.LOCAL_BUFFER)
                        fun_call.addArgument(cpu_sym_id.clone());
                    else
                        fun_call.addArgument(new UnaryExpression(UnaryOperator.ADDRESS_OF, cpu_sym_id.clone()));
                    Expression size_k = cpu_sym_id.clone();
                    size_k.setParens(true);
                    fun_call.addArgument(createDecls.createSizeofExpression(size_k));
                }
                
                else {
                    fun_call.addArgument(size.clone());
                }
                
            }

            if((from == DpuCopy.DPU_BROADCAST_TO) || (from == DpuCopy.DPU_BROADCAST_FROM)) {
                fun_call.addArgument(new NameID("DPU_XFER_DEFAULT"));
            }
            
            return fun_call;
            
        }
        
        
        ForLoop createForLoop(String it_name, Statement for_body) {
            IDExpression iter_id=new NameID(it_name);
            
            ForLoop for_loop = null;
            
                                                                
            VariableDeclarator declr = new VariableDeclarator(iter_id);
            //Initializer init = new Initializer(new IntegerLiteral(0));
            Initializer init = new Initializer(new IntegerLiteral(0));
            
            declr.setInitializer(init);
            VariableDeclaration v_decl = new VariableDeclaration(Specifier.LONG, declr);
            DeclarationStatement init_iter = new DeclarationStatement(v_decl);
            
            IDExpression limit_id = new NameID("Work_Group_Id_Limit");
            IDExpression jump_limit = new NameID("ss_nr_dpus");
            
            Expression condition = new BinaryExpression(iter_id.clone(), BinaryOperator.COMPARE_LT, limit_id);
            Expression step = new BinaryExpression(iter_id.clone(), AssignmentOperator.ADD, jump_limit);
//             Expression step = new UnaryExpression(UnaryOperator.POST_INCREMENT, iter_id.clone());
            //Statement body = for_body.clone();
            for_loop = new ForLoop(init_iter, condition, step, for_body);
            
            return for_loop;
            
        }
        
        static IDExpression getIteratorVariable(PolyStart start, int dim, String pad) {
            IDExpression iter_id = null;
            if(dim>8) {
                System.out.println("ERROR: dimension provided is a value not supported");
            }
            
//             String[] id_names = new String[]{"ss_tx", "ss_ty", "ss_tz", "ss_gx", "ss_gy", "ss_gz", "ss_bx", "ss_by", "ss_bz"};
            String[] id_names = new String[]{"tx", "ty", "tz", "gx", "gy", "gz", "bx", "by", "bz"};
            String pre = "";
            if(pad=="")
                iter_id = new NameID(pre+id_names[(int)(dim+start.getValue())]);
            else
                iter_id = new NameID(pre+id_names[(int)(dim+start.getValue())]+"_"+pad);
            
            return iter_id;
        }
        
        static IDExpression getIteratorVariable(PolyStart start, int dim, String pre, String post) {
            IDExpression iter_id = null;
            if(dim>8) {
                System.out.println("ERROR: dimension provided is a value not supported");
            }
            
//             String[] id_names = new String[]{"ss_tx", "ss_ty", "ss_tz", "ss_gx", "ss_gy", "ss_gz", "ss_bx", "ss_by", "ss_bz"};
            String[] id_names = new String[]{"tx", "ty", "tz", "gx", "gy", "gz", "bx", "by", "bz"};
            
            if(post==null)
                iter_id = new NameID(pre+id_names[(int)(dim+start.getValue())]);
            else
                iter_id = new NameID(pre+id_names[(int)(dim+start.getValue())]+post);
            
            return iter_id;
        }
        
        Expression getIteratorLimit(WorkHorseArray wha, PolyStart start, int dim) {
//             String[] limit_names = new String[]{"ss_tx_end", "ss_ty_end", "ss_tz_end", "ss_gx_end", "ss_gy_end", "ss_gz_end"};
            
            switch(start) {
                case T_START:
                    return wha.get_local_limit(dim);
                    
                case G_START:
                    return wha.get_global_limit(dim);
                    
                case B_START: {
//                     System.out.println(wha.get_global_limit(dim).toString());
//                     System.out.println(wha.get_local_limit(dim).toString());
                    return new BinaryExpression(wha.get_global_limit(dim), BinaryOperator.DIVIDE,  wha.get_local_limit(dim));
                    
                }
            }

            return null;
            
        }
        
        
        
        
        
        List<Statement> createBlockLimits(WorkHorseArray wha, int[] order, boolean create_init) {
            int len = order.length;
            List<Statement> limit_stmts = new ArrayList<Statement>();
            Expression b_expr = null;
            Expression t_expr = null;
            Expression g_expr = null;

            for(int i=0; i<len; i++) {
                IDExpression b_left = getIteratorVariable(PolyStart.B_START, order[i], "limit");
                
//                 IDExpression b_left_r = getIteratorVariable(PolyStart.B_START, order[i], "", "l");
                IDExpression t_left_r = getIteratorVariable(PolyStart.T_START, order[i], "", "l");
//                 IDExpression g_left_r = getIteratorVariable(PolyStart.G_START, order[i], "", "l");
                
                Expression t_limit = getIteratorLimit(wha, PolyStart.T_START, order[i]);
                Expression g_limit = getIteratorLimit(wha, PolyStart.G_START, order[i]);
                // System.out.println("T_LIMIT: "+ t_limit);
                if(create_init) {
//                     limit_stmts.add(createDecls.createVariableDeclaration(b_left, null));
                    
//                     limit_stmts.add(createDecls.createArrayVariableDeclaration(b_left_r, getIteratorLimit(wha, PolyStart.B_START, order[i])));
                    limit_stmts.add(createDecls.createArrayVariableDeclaration(Specifier.LONG, t_left_r, getIteratorLimit(wha, PolyStart.T_START, order[i])));
//                     limit_stmts.add(createDecls.createArrayVariableDeclaration(g_left_r, getIteratorLimit(wha, PolyStart.G_START, order[i])));
                    
                    if(i==0) {
                        b_expr = b_left.clone();
                        g_expr = g_limit;
                        // t_expr = t_limit;
                    }
                    else {
                        b_expr = new BinaryExpression(b_expr, BinaryOperator.MULTIPLY, b_left.clone());
                        g_expr = new BinaryExpression(g_expr, BinaryOperator.MULTIPLY, g_limit);
                        // t_expr = new BinaryExpression(t_expr, BinaryOperator.MULTIPLY, t_limit);
                    }
                }
                
                Expression right = getIteratorLimit(wha, PolyStart.B_START, order[i]);
                limit_stmts.add(createDecls.createVariableDeclarationStmt(Specifier.LONG, b_left, right));
                
                
                
                //limit_stmts.add(new ExpressionStatement(new AssignmentExpression(b_left.clone(), AssignmentOperator.NORMAL, right)));
            }
            if(create_init) {
//                 AssignmentExpression a_expr = new AssignmentExpression(new NameID("Work_Group_Id_Limit"), AssignmentOperator.NORMAL, b_expr); //Num of work groups this kernel call can have
                    limit_stmts.add(createDecls.createVariableDeclarationStmt(Specifier.LONG, new NameID("Work_Group_Id_Limit"), b_expr));
                
                // t_expr = new BinaryExpression(t_expr, BinaryOperator.MULTIPLY, new NameID("Work_Group_Id_Limit"));
                    limit_stmts.add(createDecls.createVariableDeclarationStmt(Specifier.LONG, new NameID("Thread_Id_Limit"), g_expr));
                // limit_stmts.add(createVariableDeclarationStmt(new NameID("Thread_Id_Limit"), t_expr));
            }
            return limit_stmts;
            
        }
            
        //bxl[] = , byl[] = , gxl[] = , gyl[] = ,
        List<Statement> createOneBlockLimits(WorkHorseArray wha, int[] order, boolean create_init, String pre, int div_group_grid) {
            int len = order.length;
            List<Statement> limit_stmts = new ArrayList<Statement>();
            
            IDExpression b_left = null;
            if(pre=="")
                b_left = new NameID("_bd");
            else
                b_left = new NameID(pre + "bd");
            List<Expression> init_exprs = new ArrayList<Expression>();

            for(int i=0; i<len; i++) {                   
                IDExpression b_left_r = getIteratorVariable(PolyStart.B_START, order[i], "", "l");
                IDExpression g_left_r = getIteratorVariable(PolyStart.G_START, order[i], "", "l");
                
                if(create_init) {
                    Expression var_exp = getIteratorVariable(PolyStart.B_START, i, pre, "");
                    
                    if(order[i]==div_group_grid) {
                        Expression dg_end = new BinaryExpression(var_exp.clone(), BinaryOperator.ADD, OneArith.minusOne(new NameID("N_MULTI_WGS")) );
                        dg_end = new FunctionCall(new NameID("ss_min"), dg_end, OneArith.minusOne(getIteratorVariable(PolyStart.B_START, order[i], "", "_limit")));
                        limit_stmts.add(createDecls.createArrayVariableDeclaration(Specifier.LONG, b_left_r, var_exp, dg_end ));
                    }
                    else
                       limit_stmts.add(createDecls.createArrayVariableDeclaration(Specifier.LONG, b_left_r, var_exp, var_exp.clone()) ); 
                    
                    // limit_stmts.add(createDecls.createArrayVariableDeclaration(g_left_r, new BinaryExpression(wha.getIteratorLimitExpr(PolyStart.T_START, order[i]), BinaryOperator.MULTIPLY, getIteratorVariable(PolyStart.B_START, order[i], pre, "")) ));
                    
                    Expression left_ind_g = new BinaryExpression(wha.getIteratorLimitExpr(PolyStart.T_START, order[i]), BinaryOperator.MULTIPLY, wha.getIteratorStartExpr(PolyStart.B_START, order[i]));
                    Expression right_ind_g = new BinaryExpression(wha.getIteratorLimitExpr(PolyStart.T_START, order[i]), BinaryOperator.MULTIPLY, wha.getIteratorLimitExpr(PolyStart.B_START, order[i]));
                    right_ind_g = new BinaryExpression(right_ind_g, BinaryOperator.ADD, OneArith.minusOne(wha.getIteratorLimitExpr(PolyStart.T_START, order[i])));
                    limit_stmts.add(createDecls.createArrayVariableDeclaration(Specifier.LONG, g_left_r, left_ind_g, right_ind_g));
                    init_exprs.add(var_exp.clone());
//                     limit_stmts.add(createDecls.createArrayVariableDeclaration(b_left_r, var_exp, new BinaryExpression(var_exp.clone(), BinaryOperator.ADD, new IntegerLiteral(1)) ));
                }
                
                
                
                //limit_stmts.add(new ExpressionStatement(new AssignmentExpression(b_left.clone(), AssignmentOperator.NORMAL, right)));
            }
            

            if(init_exprs.size()>0) {
                limit_stmts.add(new DeclarationStatement(createDecls.createArrayVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), b_left, init_exprs)));
            }
            
            limit_stmts.add(createDecls.createVariableDeclarationStmt(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), DPU_MRAM_OFFSET.clone(), new IntegerLiteral(0)));
            limit_stmts.add(createDecls.createVariableDeclarationStmt(Specifier.INT, new NameID("ss_temp_offset"), new IntegerLiteral(0)));
            limit_stmts.add(new DeclarationStatement(createDecls.createCodeAnnotDecl("int (*ss_SE)[4];")));
            limit_stmts.add(createDecls.createVariableDeclarationStmt(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG)), new NameID("ss_temp_length"), new IntegerLiteral(0)));


            return limit_stmts;
        }
        
        List<Statement> createIncBlockLimits(WorkHorseArray wha, int[] order, String pre, int div_group_grid) {
            int len = order.length;
            List<Statement> limit_stmts = new ArrayList<Statement>();
            Expression expr = null;
            Expression prev_right = null;
            IDExpression prev_left = null;

            IDExpression diff = getIteratorVariable(PolyStart.B_START, div_group_grid, "", "l");
            Declaration n_multi_wgs_decl = createDecls.createCodeAnnotDecl("int n_multi_wgs = "+diff.getName()+"[1]-"+diff.getName()+"[0]+1");

            limit_stmts.add(new DeclarationStatement(n_multi_wgs_decl));
            
            for(int i=0; i<len; i++) {
                    
                IDExpression left = getIteratorVariable(PolyStart.B_START, order[i], pre, null);

                Expression left_inc = null;
                if(i==div_group_grid)
                    left_inc = new BinaryExpression(left, BinaryOperator.ADD, new NameID("n_multi_wgs"));
                else
                    left_inc = new BinaryExpression(left, BinaryOperator.ADD, new IntegerLiteral(1));
                //Expression right = getIteratorLimit(wha, PolyStart.B_START, order[i]);
                Expression right = getIteratorVariable(PolyStart.B_START, order[i], "limit"); 
                
                if(i==0) {
                    expr = left.clone();
                    expr.setParens(false);
                }
                else {
                    Expression add_expr = new BinaryExpression(left.clone(), BinaryOperator.MULTIPLY, prev_right);
                    expr = new BinaryExpression(add_expr,BinaryOperator.ADD, expr);
                    expr.setParens(false);
                }
                prev_right = right;
                right = new BinaryExpression(left_inc, BinaryOperator.MODULUS, right.clone());
                Statement expr_stmt = new ExpressionStatement(new AssignmentExpression(left.clone(), AssignmentOperator.NORMAL, right));
                if(i==0)
                    limit_stmts.add(expr_stmt);
                else {
                    Expression condition = new BinaryExpression(prev_left.clone(), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0));
                    limit_stmts.add(new IfStatement(condition, expr_stmt));         // checking if bx==0 then increment by like that 
                }
                prev_left = left;
                
            }
            Expression left = new NameID(pre+"Work_Group_Id");
            //limit_stmts.add(new ExpressionStatement(new AssignmentExpression(left, AssignmentOperator.NORMAL, expr)));
            limit_stmts.add(new ExpressionStatement(new AssignmentExpression(left, AssignmentOperator.ADD, new NameID("n_multi_wgs"))));
            
            return limit_stmts;
        }
        
        Statement createifLimitCondition(String pre) {
            Expression left = new NameID(pre+"Work_Group_Id");
            //limit_stmts.add(new ExpressionStatement(new AssignmentExpression(left, AssignmentOperator.NORMAL, expr)));
            //limit_stmts.add(new ExpressionStatement(new AssignmentExpression(left, AssignmentOperator.ADD, new NameID("N_MULTI_WGS"))));
            Expression cond = new BinaryExpression(new NameID(pre+"Work_Group_Id"), BinaryOperator.COMPARE_LT, new NameID("Work_Group_Id_Limit"));
            
            CompoundStatement if_body = new CompoundStatement();
            if(pre != "")
                if_body.addStatement(new BreakStatement());
            else {
                IDExpression copy_id = new NameID("INPUT");
                // if_body.addStatement(new ExpressionStatement(new AssignmentExpression(copy_id.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(-1)) ));
                if_body.addDeclaration(createDecls.createVariableDeclaration(Specifier.LONG, copy_id, new IntegerLiteral(-1)));
                if_body.addStatement(makeDpuAssert( createDpuCopy(copy_id.clone(), new IntegerLiteral(0), copy_id.clone(), null, null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO) ));
                if_body.addStatement(new ContinueStatement());
            }
                
            IfStatement if_stmt = new IfStatement(new UnaryExpression(UnaryOperator.LOGICAL_NEGATION, cond.clone()), if_body);
            return if_stmt;
        }
        WhileLoop createWhileBlockLoop(WorkHorseArray wha, int[] order, CompoundStatement dpu_body, boolean create_init, int div_group_grid) {
            int len = order.length;
            
            Expression expr = null;
            
            
            List<Statement> inner_stmts = createIncBlockLimits(wha, order, "", div_group_grid);
            List<Statement> one_block_inner_stmts = createOneBlockLimits(wha, order, create_init, "", div_group_grid);
            
            // CompoundStatement block = new CompoundStatement();
            
            CompoundStatement while_body = new CompoundStatement();
            while_body.addStatement(createDpuForEach());
            while_body.addStatement(dpu_body);
            
//             dpu_body.addStatement(createifLimitCondition(""));
            DFIterator df_iter_body = new DFIterator(dpu_body);
            boolean visit_comp_body = false;
            Statement first_stmt = null;
            
            while(df_iter_body.hasNext()) {
                Object obj = df_iter_body.next();
                if((!visit_comp_body) && (obj instanceof CompoundStatement))
                {
                    visit_comp_body = true;
                    continue;
                }
                if(obj instanceof Statement) {
                    first_stmt = (Statement)obj;
                    break;
                }
            }
            
            dpu_body.addStatementBefore(first_stmt, createifLimitCondition(""));
            for(Statement st: one_block_inner_stmts) {
                if(st instanceof DeclarationStatement) {
                    dpu_body.addDeclaration(((DeclarationStatement)st).getDeclaration().clone());
                }
                else
                    dpu_body.addStatement(st);
            }

            String code_wf_range = "printf(\"WG id:%ld bxl: %ld %ld gxl: %ld %ld\\n\", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);";
            Statement print_wf_range= createCodeAnnotStmt(code_wf_range);
            dpu_body.addStatement(print_wf_range);
            
            
            
            
            //dpu_body.addStatement(new ExpressionStatement(createDpuCopy(new NameID("b_dim"), new IntegerLiteral(0), new NameID("b_dim"), null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_TO)));
            
            for(Statement st: inner_stmts) {
                if(st instanceof DeclarationStatement) {
                    dpu_body.addDeclaration(((DeclarationStatement)st).getDeclaration().clone());
                }
                else
                    dpu_body.addStatement(st);
            }
            
            
            
            FunctionCall func_call = new FunctionCall(new NameID("dpu_launch"));
            func_call.addArgument(new NameID(dpu_set_name));
            func_call.addArgument(new NameID("DPU_SYNCHRONOUS"));
            
            while_body.addStatement(new ExpressionStatement(func_call));
            Expression cond = new BinaryExpression(new NameID("Work_Group_Id"), BinaryOperator.COMPARE_LT, new NameID("Work_Group_Id_Limit"));
//             IfStatement if_stmt = new IfStatement(new UnaryExpression(UnaryOperator.LOGICAL_NEGATION, cond.clone()), new BreakStatement());
//             dpu_body.addStatement(if_stmt);
            
            WhileLoop wh = new WhileLoop(cond, while_body);
            // prev_while_dpu_block = while_body;
            
            // block.addStatement(wh);
            // return block;
            return wh;
            
        }

        Statement createBlock(Statement stmt) {
            CompoundStatement block = new CompoundStatement();
            block.addStatement(stmt);
            return block;
        }
	
        
        /**/
        
        Statement createAllocDpu() {
//             List<Statement> alloc_stmts = new ArrayList<Statement>();
            String code = "\nstruct dpu_set_t ss_dpu_set, ss_dpu;\n"
                            +"char* ss_profile = \"backend=simulator\";\n"
                            +"DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));\n"
                            +"uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; \n"
                            +"dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); \n"
                            +"printf(\"DPU_ALLOCATED = %d\\n\", ss_nr_dpus); \n"
                            +"dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); \n"
                            +"printf(\"DPU_RANKS = %d\\n\", ss_nr_ranks); \n\n";
            
            return createDecls.createCodeAnnotStmt(code);
        }
        
        Statement createFreeDpu() {
            String code = "DPU_ASSERT( dpu_free(ss_dpu_set) );\n";
            
            return createDecls.createCodeAnnotStmt(code);
        }
        
        Statement createCountCyclesDpu(int k_ind) {
            String code = "uint64_t dpu_id = 0;\n"+
                        "int64_t cycles[ss_nr_dpus];\n"+
                        "int64_t barrier_count[ss_nr_dpus][N_TASKLETS];\n"+
                        "int64_t work_done[ss_nr_dpus][N_TASKLETS];\n"+
                        "//int64_t T_TD_START[ss_nr_dpus][N_TASKLETS];\n"+
                        "//int64_t T_TD_END[ss_nr_dpus][N_TASKLETS];\n"+
                        "//int64_t N_WG_ID[ss_nr_dpus][N_TASKLETS];\n"+
                        //"int64_t b_dim[ss_nr_dpus][3];\n"+
                        "int ss_used_ndpus = (Work_Group_Id_Limit+N_MULTI_WGS-1)/N_MULTI_WGS;\n"+
                        "int ss_ndpus = (ss_nr_dpus < ss_used_ndpus) ? ss_nr_dpus : ss_used_ndpus;\n\n"+
                        "for(int ss_c=0; ss_c<ss_nr_dpus; ss_c++)\n"+
                        "\tcycles[ss_c]=0;\n"+

                        "int ss_iter_dpu;\n"+
                        "DPU_FOREACH(ss_dpu_set, ss_dpu, ss_iter_dpu) {\n"+
                        "if(ss_iter_dpu>ss_ndpus)\n"+
                        "break;\n"+
                        "//DPU_ASSERT(dpu_log_read(ss_dpu, stdout));\n"+
                        "DPU_ASSERT(dpu_copy_from(ss_dpu, \"cycles\", 0, &cycles[dpu_id], sizeof(int64_t)));\n"+
                        "//DPU_ASSERT(dpu_copy_from(ss_dpu, \"barrier_count\", 0, &barrier_count[dpu_id], sizeof(barrier_count[0])));\n"+
                        "DPU_ASSERT(dpu_copy_from(ss_dpu, \"WORK_DONE\", 0, &work_done[dpu_id], sizeof(work_done[0])));\n"+
                        "//DPU_ASSERT(dpu_copy_from(ss_dpu, \"T_TD_START\", 0, &T_TD_START[dpu_id], sizeof(T_TD_START[0])));\n"+
                        "//DPU_ASSERT(dpu_copy_from(ss_dpu, \"T_TD_END\", 0, &T_TD_END[dpu_id], sizeof(T_TD_END[0])));\n"+
                        "//DPU_ASSERT(dpu_copy_from(ss_dpu, \"N_WG_ID\", 0, &N_WG_ID[dpu_id], sizeof(N_WG_ID[0])));\n"+
                        //"DPU_ASSERT(dpu_copy_from(ss_dpu, \"_bd\", 0, &b_dim[dpu_id], sizeof(b_dim[0])));\n"+
                        "dpu_id++;\n"+
                        "}\n\n"+
                        "uint64_t total_cycles =0;\n"+
                        "int32_t max_count = 0;\n\n"+
                        
                        "for(int i=0; i<ss_ndpus; i++) {\n"+
                        // "for(int i=0; i<ss_nr_dpus; i++) {\n"+
                        "printf(\"%ld \", cycles[i]);\n"+
                        "total_cycles+=cycles[i];\n"+
                        "if(max_count<cycles[i])\n"+
                        "\tmax_count = cycles[i];\n"+
                        "}\n\n"+

                        "printf(\"Work_Group_Id: %ld\\n\", Work_Group_Id);\n"+
                        "printf(\"kernel - Total cycles = %\" PRId64 \"\\n\", total_cycles);\n"+
                        "printf(\"kernel - Max cycles = %\" PRId32 \"\\n\", max_count);\n\n"+
                        "printf(\"----------\\n\");\n"+
                        "int trace_values[5] = {N_TASKLETS, N_MULTI_WGS, PARTITION_DIM_GRID, PARTITION_DIM_WG, max_count};\n"+
                        "note_down(trace_values, "+k_ind+");\n\n"+
                        "for(int i=0; i<ss_ndpus; i++) {\n"+
                        // "printf(\"Dpu %d:\\t\",i);\n"+
                        "for(int j=0; j<N_TASKLETS; j++) {\n"+
                        "//printf(\"%d -> T:%ld-%ld, WG:%ld, B:%ld, WD:%ld \", j, T_TD_START[i][j], T_TD_END[i][j], N_WG_ID[i][j], barrier_count[i][j], work_done[i][j]);\n"+
                            "if(work_done[i][j] == 0) \n"+
                            "\t\tprintf(\"[ERROR] WORK NOT DONE : Dpu %d:\\n\", i);\n"+
                        "}\n"+
                        "//printf(\"\\n\");\n"+
                        "}\n";
                        
            return createCodeAnnotStmt(code);
        }
        
        Statement createCodeAnnotStmt(String code) {
            CodeAnnotation c_annot = new CodeAnnotation(code);
            AnnotationStatement annot_stmt = new AnnotationStatement(c_annot);
            return annot_stmt;
        }
        
        Statement createDpuLoad(Expression kernel_id) {
        
            FunctionCall func_call = new FunctionCall(new NameID("dpu_load"));
            
            func_call.addArgument(new NameID(dpu_set_name));
            func_call.addArgument(kernel_id.clone());
            func_call.addArgument(new NameID("NULL"));
            return makeDpuAssert(func_call);
            
        }
        
        Statement createDpuLaunch() {
            FunctionCall func_call = new FunctionCall(new NameID("dpu_launch"));
            func_call.addArgument(new NameID(dpu_set_name));
            func_call.addArgument(new NameID("DPU_SYNCHRONOUS"));
            return makeDpuAssert(func_call);
        }
        
        Statement createDpuForEach() {
            return createCodeAnnotStmt(createDpuLoop().toString());
        }
        
        Expression createDpuLoop() {
            FunctionCall fun_call = new FunctionCall(new NameID("DPU_FOREACH"));
            fun_call.addArgument(new NameID(dpu_set_name));
            fun_call.addArgument(new NameID(dpu_name));
            return fun_call;
        }
        
        void convertReadBuffer(FunctionCall func_call, KernelRepr _kernel) {
            // System.out.println("ConvertReadBuffer");
            Expression gpu_buf = func_call.getArgument(1);
            Expression cpu_buf = func_call.getArgument(5);
            Expression offset = func_call.getArgument(3);
            Expression size_bytes = func_call.getArgument(4);
            
            System.out.println("GPU BUF: "+ gpu_buf); 
            if(_kernel == null) {
                System.out.println("Error: Kernel not created");
                return;
            }  
            int arg_index = _kernel.arguments_map.get((IDExpression)gpu_buf);
            TranslationUnitKernelInfo t_info = translation_unit_info.get(_kernel.get_kernel_name_string());
            CompoundStatement cur_read_comp_block = _kernel.read_comp_block;
            Statement cur_read_comp_block_ref = _kernel.read_comp_block_ref;
//             int arg_index = t_info.primitive.get_parameter_id((IDExpression)gpu_buf);
            
            List< Triple<Statement, Integer, AccessType> > offset_start_stmts = _kernel.get_offset_expression(arg_index, gpu_buf, cur_read_comp_block, "p_", true, true); // numder of index types that are gettig used to access this argument
            List< Triple<Statement, Integer, AccessType> > offset_end_stmts = _kernel.get_offset_expression(arg_index, gpu_buf, cur_read_comp_block, "p_", false, true);
            
            
            if((offset_start_stmts == null) || (offset_start_stmts.size()==0))
                    return;

            CompoundStatement dpu_copy_from_stmts = new CompoundStatement();
            // CreateCopyReadStatements(t_info, arg_index, (IDExpression)gpu_buf, _kernel, cpu_buf, offset_start_stmts, offset_end_stmts, DpuCopy.DPU_COPY_FROM, dpu_copy_from_stmts);
            CreateCopyReadStatements(t_info, _kernel, arg_index, (IDExpression)gpu_buf, cpu_buf, DpuCopy.DPU_COPY_FROM, dpu_copy_from_stmts);

            cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, dpu_copy_from_stmts);

            

            // int num_off_exprs = offset_start_stmts.size();


            // Iterator< Triple<Statement, Integer, AccessType> > offset_start_iterator = offset_start_stmts.iterator();
            // Iterator< Triple<Statement, Integer, AccessType> > offset_end_iterator = offset_end_stmts.iterator(); 
            
            // List<VariableDeclarator> var_declr = new ArrayList<VariableDeclarator>();

            // if(num_off_exprs==0)
                // return;

            



            // Expression mram_offset_init = new AssignmentExpression(DPU_MRAM_OFFSET.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0));

            // cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, new ExpressionStatement(mram_offset_init));

            // // func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), new IntegerLiteral(0), DPU_MRAM_OFFSET.clone(), null, ArgType.NORMAL, null, DpuCopy.DPU_COPY_FROM);
            // // cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, makeDpuAssert(func_call));

            // Expression symbol_offset = DPU_MRAM_OFFSET.clone();
            // AssignmentOperator assign_op = AssignmentOperator.NORMAL;
            // Expression size_id = new NameID("p_"+((IDExpression)gpu_buf).getName()+"_size");
            
            // System.out.println(_kernel.get_kernel_name_string() + " (" + arg_index + ") : Number of READ offset expressions: " + num_off_exprs);
            
            // IDExpression ele_size_id = new NameID("p_"+((IDExpression)gpu_buf).getName()+"_ele_size");
            
            // Statement ele_size = createElementSizeStmt(cpu_buf, ele_size_id);
            // cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, ele_size);
            
            // for(int i=0; i<num_off_exprs; i++) {

            //     Triple<Statement, Integer, AccessType> pair_start = offset_start_iterator.next();
            //     ExpressionStatement offset_start_expr = (ExpressionStatement)pair_start.getFirst();
            //     cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, offset_start_expr);
            //     AssignmentExpression start_expr = (AssignmentExpression)offset_start_expr.getExpression();

            //     Triple<Statement, Integer, AccessType> pair_end = offset_end_iterator.next();
            //     ExpressionStatement offset_end_expr = (ExpressionStatement)pair_end.getFirst();
            //     AssignmentExpression end_expr = (AssignmentExpression)offset_end_expr.getExpression();
                
            //     if(start_expr.getRHS().compareTo(end_expr.getRHS())==0) {
            //         end_expr.getRHS().swapWith( new BinaryExpression(end_expr.getRHS().clone(), BinaryOperator.ADD, new IntegerLiteral(2)) );
            //     }

            //     cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, offset_end_expr);

            //     Expression size_expr = _kernel.get_offset_expression_size(gpu_buf, pair_start.getSecond(), "p_");

            //     if(i!=0) {
            //         symbol_offset = new BinaryExpression(DPU_MRAM_OFFSET.clone(), BinaryOperator.ADD, new BinaryExpression(size_id, BinaryOperator.MULTIPLY, ele_size_id));
            //         assign_op = AssignmentOperator.ADD;
            //     }
                
                
            //     func_call = createDpuCopy(t_info.primitive.get_parameter_idexp(arg_index), symbol_offset, cpu_buf, size_expr, ArgType.GLOBAL_BUFFER, gpu_buf,  DpuCopy.DPU_COPY_FROM);
            //     cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, makeDpuAssert(func_call));
                
            //     cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, new ExpressionStatement(new AssignmentExpression(size_id.clone(), assign_op, size_expr.clone())));

            //     if(i==num_off_exprs-1) {
            //         Expression dpu_mram_offset_modify = new AssignmentExpression(DPU_MRAM_OFFSET.clone(), AssignmentOperator.ADD, new BinaryExpression(size_id.clone(), BinaryOperator.MULTIPLY, ele_size_id.clone()));
            //         cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, new ExpressionStatement(dpu_mram_offset_modify));
            //     }

            // }
            // var_declr.add(new VariableDeclarator(new NameID("p_"+((IDExpression)gpu_buf).getName()+"_offset")) );
            // var_declr.add(new VariableDeclarator(new NameID("p_"+((IDExpression)gpu_buf).getName()+"_size")) );
            // var_declr.add(new VariableDeclarator(new NameID("p_"+((IDExpression)gpu_buf).getName()+"_end")) );
            // var_declr.add(new VariableDeclarator(ele_size_id.clone()) );
            // cur_read_comp_block.addStatementBefore(cur_read_comp_block_ref, createCodeAnnotStmt("\n"));
            
            // cur_read_comp_block.addDeclaration(new VariableDeclaration(Specifier.INT, var_declr));
            
        }
        
        boolean swap_expressions_IR()
        {
            // Iterator<Expression> ref_exprs = ref_swap_exprs.iterator();

            for(Expression ref_expr: ref_swap_exprs) {
                ref_expr.swapWith(CL_SUCCESS.clone());
            }

            return true;
        }
        
        void swap_expression(Expression current_expr, Expression new_expr) 
        {
            //changed_stmt = if_stmt;
            if(new_expr==null) {
                ref_swap_exprs.add(current_expr);
            }
            
            // CHECK: else part : FUTURE
        }

        boolean swap_statements_IR()
        {
            Iterator<Statement> ref_stmts = ref_swap_stmts.iterator();
            Iterator<Statement> comp_swap_stmts_it = comp_swap_stmts.iterator();
            //Iterator<Boolean> swap_to_before_it = is_swap_next.iterator();
            
            //Iterator<Object> comp_stmts_swap_add = comp_swap_stmts.iterator();
            int add_stmts = 0;

            for(Object st: add_swap_stmts) {
                ++add_stmts;
                int error=0;
                Statement new_stmt = (Statement)st;
                Statement ref_st = (Statement) ref_stmts.next();
                CompoundStatement c_st = (CompoundStatement) comp_swap_stmts_it.next();
                
                ref_st.swapWith(new_stmt);
                
                
            }
            return true;
        }
        
        boolean add_statements_IR() 
        {
            Iterator<Statement> ref_to_stmts_it = ref_to_stmts.iterator();
            Iterator<Statement> comp_add_stmts_it = comp_add_stmts.iterator();
            //Iterator<Boolean> add_to_before_it = is_add_next.iterator();
            Iterator<Boolean> is_add_before_it = is_add_before.iterator();

            for(Object st: add_to_stmts) {

                int error=0;
                Statement new_stmt = (Statement)st;
                Statement ref_st = (Statement) ref_to_stmts_it.next();
                CompoundStatement c_st = (CompoundStatement) comp_add_stmts_it.next();
                boolean is_before = is_add_before_it.next();
                try{
                if(is_before)
                    c_st.addStatementBefore(ref_st, new_stmt);
                else
                    c_st.addStatementAfter(ref_st, new_stmt); 
                }
                catch(IllegalArgumentException e) {
                    System.out.println(new_stmt.toString());
                    throw e;
                }

            }
            return true;
        }
        
        boolean rm_statements_IR() 
        {
            //Iterator<Statement> ref_to_stmts_it = rm_to_stmts.iterator();
            Iterator<Statement> comp_add_stmts_it = comp_rm_stmts.iterator();
            //Iterator<Boolean> add_to_before_it = is_add_next.iterator();

            for(Object st: rm_to_stmts) {

                int error=0;
                Statement rm_st = (Statement) st;
                CompoundStatement c_st = (CompoundStatement) comp_add_stmts_it.next();
                
                try {
//                     if(!(not_rm_stmt.contains(rm_st))) {
//                         System.out.println("REMOVE: "+ rm_st);
                        // System.out.println("REMOVE: "+ rm_st);
                        c_st.removeStatement(rm_st);
//                     }
                }
                catch(IllegalArgumentException e) {
                    System.out.println(rm_st.toString());
                    throw e;
                }

            }
            
//             for(Expression expr: rm_exprs) {
//                 List<Expression> rm_expressions = IRTools.findExpressions(tran_unit, expr);
//                 for(Expression sub_expr: rm_expressions) {
//                     Statement st = sub_expr.getStatement();
//                     if(st!=null && !(not_rm_stmt.contains(st))) {
//                         System.out.println("REMOVE: " +expr + " :: " + st);
// //                         st.swapWith(new NullStatement());
//                     }
// //                         st.swapWith(createCodeAnnotStmt(""));
// //                         st.swapWith(new NullStatement());
//                 }
//             }
            return true;
        }
    
        void swap_with_statement(Statement current_stmt, Statement changed_stmt, Statement current_comp_stmt) 
        {
            ref_swap_stmts.add(current_stmt);
            //changed_stmt = if_stmt;
            add_swap_stmts.add(changed_stmt);
            comp_swap_stmts.add(current_comp_stmt);
        }
        
        void add_to_statement(Statement current_stmt, Statement changed_stmt, Statement current_comp_stmt, boolean is_before) 
        {
            //changed_stmt = if_stmt;
            ref_to_stmts.add(current_stmt);
            add_to_stmts.add(changed_stmt);
            comp_add_stmts.add(current_comp_stmt);
            is_add_before.add(is_before);
        }
        
        void rem_statement(Statement current_stmt, Statement current_comp_stmt) 
        {
            //changed_stmt = if_stmt;
            rm_to_stmts.add(current_stmt);
            comp_rm_stmts.add(current_comp_stmt);
        }

        
        
}

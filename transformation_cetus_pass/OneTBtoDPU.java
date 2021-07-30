package cetus.transforms;
import cetus.analysis.AnalysisPass;
import cetus.analysis.RangeAnalysis;
import cetus.analysis.DDGraph;
import cetus.analysis.DependenceVector;
import cetus.analysis.LoopTools;
import cetus.analysis.RangeDomain;
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
// import org.apache.commons.io.FilenameUtils;

//import cetus.exec.java;
class Pair<T1, T2>
{
    T1 p1;
    T2 p2;
    Pair(T1 p1, T2 p2) {
        this.p1 = p1;
        this.p2 = p2;
    }
    
    void setFirst(T1 p1) {
        this.p1 = p1;
    }
    void setSecond(T2 p2) {
        this.p2 = p2;
    }
    
    T1 getFirst() {
        return p1;
    }
    
    T2 getSecond() {
        return p2;
    }

    void print() {
        System.out.println("0: " + p1+ " 1: " + p2);
    }
    
}

class Triple<T1, T2, T3>
{
    T1 p1;
    T2 p2;
    T3 p3;
    Triple(T1 p1, T2 p2, T3 p3) {
        this.p1 = p1;
        this.p2 = p2;
        this.p3 = p3;
    }
    
    void setFirst(T1 p1) {
        this.p1 = p1;
    }
    void setSecond(T2 p2) {
        this.p2 = p2;
    }
    void setThird(T3 p3) {
        this.p3 = p3;
    }
    
    T1 getFirst() {
        return p1;
    }
    
    T2 getSecond() {
        return p2;
    }
    T3 getThird() {
        return p3;
    }
    
}

class Quadruple<T1, T2, T3, T4>
{
    T1 p1;
    T2 p2;
    T3 p3;
    T4 p4;
    Quadruple(T1 p1, T2 p2, T3 p3, T4 p4) {
        this.p1 = p1;
        this.p2 = p2;
        this.p3 = p3;
        this.p4 = p4;
    }
    
    void setFirst(T1 p1) {
        this.p1 = p1;
    }
    void setSecond(T2 p2) {
        this.p2 = p2;
    }
    void setThird(T3 p3) {
        this.p3 = p3;
    }
    void setFourth(T4 p4) {
        this.p4 = p4;
    }
    
    T1 getFirst() {
        return p1;
    }
    T2 getSecond() {
        return p2;
    }
    T3 getThird() {
        return p3;
    }
    T4 getFourth() {
        return p4;
    }

}

enum PolyStart {
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
}

/*class Coefficients {
    int total_coeff = 9;
    List<Expression>[] co_eff = new ArrayList[total_coeff];
    long set_coeff = 0;
    
    Coefficients() {
        for(int i=0; i<total_coeff; i++) {
            co_eff[i] = new ArrayList<Expression>();
        }
    }
    
    
        
//     void set_coefficient(Expression expr, PolyStart start, int dim) {
//         co_eff[start.getValue()+dim].add(expr);
//         set_coeff |= (1<< (start.getValue()+dim));
//     }
//     
//     List<Expression> get_coefficient(PolyStart start, int dim) {
//         return co_eff[start.getValue()+dim];
//     }
    void set_coefficient(Expression expr, int dim) {
        co_eff[dim].add(expr);
        set_coeff |= (1<< dim);
    }
    
    List<Expression> get_coefficient(int dim) {
        return co_eff[dim];
    }
    
    boolean set_coefficients(Expression expr, Primitive primitive) {    //CAUTION: pass clone of the expression.
        if(expr instanceof BinaryExpression) {
            BinaryExpression b_expr = (BinaryExpression) expr;
            if(b_expr.getOperator().equals(BinaryOperator.ADD)) {
                set_coefficients(b_expr.getLHS(), primitive);
                set_coefficients(b_expr.getRHS(), primitive);
            }
            
        }
        else {
            DFIterator<IDExpression> expr_iterator = new DFIterator<IDExpression>(expr, IDExpression.class);
            
            while(expr_iterator.hasNext()) {
                Expression sub_expr = expr_iterator.next();
                
                int res = primitive.is_primitive(sub_expr);
                if(res==-1)
                    return false;
                if(res<64) {
                    IDExpression id_iter = OneTBtoDPU.getIteratorVariable(res);
                    
                    sub_expr.swapWith(id_iter);
                }
                
            }
        }
        
        return true;
    }
}*/

class OneArith {
    static BinaryExpression minusOne(Expression expr) {
        return new BinaryExpression(expr, BinaryOperator.SUBTRACT, new IntegerLiteral(1));
    }

    static BinaryExpression addOne(Expression expr) {
        return new BinaryExpression(expr, BinaryOperator.ADD, new IntegerLiteral(1));
    }
}

class AnalyseArrayExpression {
    static Expression removeArrayAccess(Expression expr) {
        if(expr instanceof ArrayAccess) {
            return ((ArrayAccess)expr).getArrayName();
        }

        // Expression dup_expr = expr.clone();
        DFIterator<ArrayAccess> expr_iterator = new DFIterator<ArrayAccess>(expr, ArrayAccess.class);

        while(expr_iterator.hasNext()) {
            ArrayAccess acc = expr_iterator.next();
            acc.swapWith(acc.getArrayName());
        }
        return expr;
    }

    static Triple<Integer, Integer, Boolean> analyseIndex(Expression expr, Strides coeff, int cur_ind, Primitive primitive, KernelRepr _kernel) {     //remove the constant terms
        if(expr == null)
            return new Triple(0, cur_ind, true);
        if(expr instanceof Literal) {
            return new Triple(0, cur_ind, true);
        }
        else if(expr instanceof IDExpression) {
            IDExpression id_expr = (IDExpression)expr;

            if(_kernel.get_argument_index(id_expr)!=-1) {
                id_expr = primitive.total_parameter_ind_map.get(_kernel.get_argument_index(id_expr));
                 System.out.println("trans: "+id_expr + " ");
            }
            else if(id_expr.getName().startsWith("host_")) {
                id_expr = new NameID(id_expr.getName().substring(5));
                System.out.println(id_expr + " "+ primitive.is_primitive(id_expr));
            }

            if(primitive.scalars.contains(id_expr)) {
                // System.out.println("scalar: "+ expr);
                return new Triple(0, cur_ind, true);
            }
            
            int prim_ind = primitive.is_primitive(id_expr);
            
            if(prim_ind>=0) {

                if(prim_ind<64 || prim_ind>=96) {
                // if(cur_ind!=0)
                    // return new Pair(-1, cur_ind);

                    int index = coeff.set_stride((IDExpression)expr, new IntegerLiteral(1));    //CAUTION: dont change it to id_expr
                    // System.out.println("Index of "+ expr + "= "+ index);
                    expr.swapWith(new IntegerLiteral(1));

                    return new Triple(1, index, false);
                }

                // else if(prim_ind>=96) {
                //     int index = coeff.set_stride((IDExpression)expr, new IntegerLiteral(1));    //CAUTION: dont change it to id_expr
                //     expr.swapWith(new IntegerLiteral(1));

                //     return new Triple(1, index, false);
                // }
            }

        }
        else if(expr instanceof UnaryExpression) {
            UnaryOperator u_op = ((UnaryExpression)expr).getOperator();

            if((u_op == UnaryOperator.MINUS)        //check in detail- will minus work? <- DOUBT
                || (u_op == UnaryOperator.PLUS))
                return analyseIndex(((UnaryExpression)expr).getExpression(), coeff, cur_ind, primitive, _kernel);

            return new Triple(-1, cur_ind, false);
        }
        if(!(expr instanceof BinaryExpression))
            return new Triple(-1, cur_ind, false);

        BinaryExpression bin_expr = (BinaryExpression)expr;

        BinaryOperator bin_op = bin_expr.getOperator();
        Expression bin_lhs = bin_expr.getLHS();
        Expression bin_rhs = bin_expr.getRHS();

        Triple<Integer, Integer, Boolean> _lhs = analyseIndex(bin_lhs, coeff, cur_ind, primitive, _kernel);
        Triple<Integer, Integer, Boolean> _rhs = analyseIndex(bin_rhs, coeff, cur_ind, primitive, _kernel);

        if(_lhs.getFirst()==-1 || _rhs.getFirst()==-1 )
            return new Triple(-1, cur_ind, false);

        else if( bin_op== BinaryOperator.ADD || bin_op==BinaryOperator.SUBTRACT) {

            if(_lhs.getFirst()==1) {
                // change_stride(_lhs.getSecond(), bin_lhs);
                cur_ind = _lhs.getSecond();
                if(cur_ind!=0 && _lhs.getThird()) {
                    // System.out.print("checkl: "+ bin_lhs);
                    Expression swap_expr = IRTools.findExpression(bin_lhs, new IntegerLiteral(1));
                    if(swap_expr!=null)
                        swap_expr.swapWith(new IntegerLiteral(0));
                    // System.out.println(" "+ bin_lhs);
                }
            }

            if(_rhs.getFirst()==1) {
                // change_stride(_rhs.getSecond(), bin_rhs);
                cur_ind = _rhs.getSecond();
                // System.out.print("checkr: "+ bin_rhs);
                if(cur_ind!=0 && _rhs.getThird()) {
                    Expression swap_expr = IRTools.findExpression(bin_rhs, new IntegerLiteral(1));
                    if(swap_expr!=null)
                        swap_expr.swapWith(new IntegerLiteral(0));
                    // System.out.println(" "+ bin_rhs);
                }
                // bin_rhs.swapWith(new IntegerLiteral(0));
            }

            return new Triple((_lhs.getFirst()+_rhs.getFirst()), cur_ind, _lhs.getThird()&&_rhs.getThird());
        }
        else if(bin_op == BinaryOperator.MULTIPLY) {
            if(_lhs.getFirst()+_rhs.getFirst()>1)
                return new Triple(-1, cur_ind, false);

            if(_lhs.getFirst()==1) {
                cur_ind = _lhs.getSecond();
                coeff.append_stride(_lhs.getSecond(), bin_rhs.clone());
            }
            if(_rhs.getFirst()==1) {
                cur_ind = _rhs.getSecond();
                coeff.append_stride(_rhs.getSecond(), bin_lhs.clone());
            }

            return new Triple((_lhs.getFirst()+_rhs.getFirst()), cur_ind, _lhs.getThird()&&_rhs.getThird());
        }

        // if(_lhs+_rhs > 1)
            // return -1;
        return new Triple(-1, cur_ind, false);


        // if(_lhs==1 && _rhs==1)
            // return 1;
        // return 0;
    }
}

class Strides {
    Expression[] stride = {null, null, null};
    IDExpression[] s_id = {null, null, null};
    int over = 0;
    int is_strided = 0; // 0 - row access, 1-col access, -1-index pattern not row or col wise
    Expression offset = new IntegerLiteral(0);
    int pattern = 0;
    int[] pattern_ind = {-1, -1, -1};


    boolean check() {
        if(over<0 || over>2)
            return false;
        if(is_strided==-1)
            return false;
        return true;
    }

    boolean check(int index) {
        if(index<0 || index>=over)
            return false;
        if(is_strided==-1)
            return false;
        return true;
    }

    int set_stride(IDExpression _id, Expression expr) {
        if(!check())
            return -1;

        stride[over] = expr;
        s_id[over] = _id;
        return ++over;
    }

    int set_only_stride(Expression expr) {
        if(!check())
            return -1;

        stride[over] = expr;
        return ++over;
    }

    int set_only_stride_id(IDExpression _id) {
        if(!check())
            return -1;
        
        s_id[over] = _id;
        return ++over;
    }

    void set_offset(Expression expr) {
        if(expr!=null)
            // offset = expr;
            offset = Symbolic.simplify(expr);
    }

    boolean change_stride(int index, Expression _stride) {
        if(!check(index-1))
            return false;
        stride[index-1] = _stride;
        return true;
    }

    boolean append_stride(int index, Expression _stride) {
        if(!check(index-1))
            return false;
        stride[index-1] = new BinaryExpression(stride[index-1], BinaryOperator.MULTIPLY, _stride);
        return true;
    }

    boolean change_id(int index, IDExpression _id) {
        if(!check(index-1))
            return false;
        s_id[index-1] = _id;
        return true;
    }

    IDExpression get_id(int index){
        if(!check(index-1))
            return null;
        return s_id[index-1];
    }

    Map<IDExpression, Long> myMap = new HashMap<IDExpression, Long>() {{
        put(new NameID("txl"), 1l<<0);
        put(new NameID("tyl"), 1l<<1);
        put(new NameID("tzl"), 1l<<2);

        put(new NameID("gxl"), 1l<<3);
        put(new NameID("gyl"), 1l<<4);
        put(new NameID("gzl"), 1l<<5);

        put(new NameID("bxl"), 1l<<6);
        put(new NameID("byl"), 1l<<7);
        put(new NameID("bzl"), 1l<<8);
    }};

    Long get_id_vector() {
        long ind_var_access = 0l;

        for(int i=0; i<over; i++)
            ind_var_access |= myMap.get(s_id[i]);

        return ind_var_access;
    }

    void simplify_all() {
        // int cols=0;
        // Literal check_literal = new IntegerLiteral(1);
        for(int i=0; i<over; i++) {
            stride[i] = Symbolic.simplify(stride[i]);
            // if(!stride[i].equals(check_literal)) {
            //     cols++;
            //     pattern_ind[1] = i;
            // }
            // else
            //     pattern_ind[0] = i;
        }
        // if(cols==1)
        //     pattern = 1;
        // else if(cols > 1)
        //     pattern = -1;   //not implemented
    }

    void get_overall_stride(IDExpression div_group_grid) {
        int cols=0;
        Literal check_literal = new IntegerLiteral(1);

        for(int i=0; i<over; i++) {
            // stride[i] = Symbolic.simplify(stride[i]);
            if(myMap.containsKey(s_id[i])) {
                if((myMap.get(s_id[i])&448l)!=0 && !s_id[i].equals(div_group_grid)) {
                    System.out.println(s_id[i] + " "+ div_group_grid);
                    System.out.println("continue");
                    continue;
                }
            }
            

            if(!stride[i].equals(check_literal)) {
                cols++;
                pattern_ind[1] = i;
            }
            else
                pattern_ind[0] = i;
        }
        if(cols==0)
            System.out.println("overall_stride = "+ 1);
        else if(cols==1) {
            pattern = 1;
            System.out.println("overall_stride = "+ stride[pattern_ind[1]]);
        }
        else if(cols > 1)
            pattern = -1;   //not implemented

    }

    void print() {
        for(int i=0; i<over; i++) {
            System.out.println(s_id[i] + " stride= " + stride[i]);
        }
        System.out.println("   offset= " + offset);
        if(pattern==0)
            System.out.println("Row access pattern");
        else if(pattern == 1)
            System.out.println("Column access pattern : "+ stride[pattern_ind[1]]);
        else
            System.out.println("pattern analysis NOT IMPLEMENTED");
    }

    Pair<Expression, Expression> createStridedExpr(boolean stride_trunc) {
        if(pattern==1) {
            IDExpression iter_id = new NameID("ss_i");
            Expression _stride = stride[pattern_ind[1]].clone();
            Expression new_stride = null;
            if(stride_trunc){
                if(pattern_ind[0]==-1)
                    _stride = new IntegerLiteral(1);
                else {
                    _stride = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(1));
                    _stride = new BinaryExpression(_stride, BinaryOperator.SUBTRACT, new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(0)));
                    _stride = OneArith.addOne(_stride);
                }
                new_stride = _stride;
                
            }
            Expression col_expr = new BinaryExpression(_stride, BinaryOperator.MULTIPLY, iter_id.clone());

            if(pattern_ind[0]!=-1){
                // Expression row_expr = new BinaryExpression(stride[pattern_index[0]].clone(), BinaryOperator.MULTIPLY, s_id[pattern_ind[0]].clone());
                // Expression row_expr = s_id[pattern_ind[0]].clone();
                Expression new_row_expr = new NameID("ss_j");
                col_expr = new BinaryExpression(col_expr, BinaryOperator.ADD, new_row_expr);
            }

            if(!stride_trunc)
                col_expr = new BinaryExpression(col_expr, BinaryOperator.ADD, offset.clone());

            return new Pair(col_expr, new_stride);
        }
        return null;
        
    }


    static Expression createMallocTemp(Specifier spec, IDExpression temp_id, Expression length) {
        return createDecls.createMallocVariableAssign(spec, temp_id, length);
    }

    Expression tempBuffLength() {
        Expression length = null;
        if(pattern==1) {
            ArrayAccess for_start = new ArrayAccess(s_id[pattern_ind[1]].clone(), new IntegerLiteral(0));
            ArrayAccess for_end = new ArrayAccess(s_id[pattern_ind[1]].clone(), new IntegerLiteral(1));
            length = OneArith.addOne(new BinaryExpression(for_end, BinaryOperator.SUBTRACT, for_start));
            
            if(pattern_ind[0]!=-1) {
                for_start = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(0));
                for_end = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(1));
                Expression length_2 = OneArith.addOne(new BinaryExpression(for_end, BinaryOperator.SUBTRACT, for_start));
                length = new BinaryExpression(length, BinaryOperator.MULTIPLY, length_2);
            }
        }
        return length;
    }

    Pair<CompoundStatement, Expression> generateColAccCopy(Specifier spec, int arg_ind, Expression arr, IDExpression gpu_arr, IDExpression temp, DpuCopy copy_direction, FunctionCall dpu_copy_call, WorkFuncRangeMap wf_range_map) {
        CompoundStatement c_st = new CompoundStatement();
        Expression new_stride = null;
        if(pattern==1) {
            System.out.println("Arg Type: "+spec);
            ForLoop for_loop = null;
            ArrayAccess for_start = new ArrayAccess(s_id[pattern_ind[1]].clone(), new IntegerLiteral(0));
            ArrayAccess for_end = new ArrayAccess(s_id[pattern_ind[1]].clone(), new IntegerLiteral(1));
            
            
            IDExpression temp_length_id = new NameID(temp.getName()+"_length");
            VariableDeclaration v_decl = createDecls.createVariableDeclaration(Specifier.INT, temp_length_id, tempBuffLength());
            c_st.addDeclaration(v_decl);



            String pre = "";
            if(copy_direction==DpuCopy.DPU_COPY_FROM || copy_direction==DpuCopy.DPU_BROADCAST_FROM) {
                pre = "p_";    
            }
            Expression _end = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_end"), new IntegerLiteral(arg_ind));
            Expression _start = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_offset"), new IntegerLiteral(arg_ind));
            // Expression _stride = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_stride"), new IntegerLiteral(arg_ind));

            Expression length_row = OneArith.addOne(new BinaryExpression(_end.clone(), BinaryOperator.SUBTRACT, _start.clone()));
            Expression cond = new BinaryExpression(length_row, BinaryOperator.COMPARE_GT, temp_length_id.clone());

            CompoundStatement if_c_st = new CompoundStatement();
            Expression modify_length = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, KernelRepr.make8ByteAligned(temp_length_id.clone(), createDecls.createSizeofExpression(new NameID(spec.toString()))));
            if_c_st.addStatement(new ExpressionStatement(modify_length));
            if_c_st.addStatement(new ExpressionStatement(createMallocTemp(spec, temp, temp_length_id.clone())) );
            

            AssignmentExpression change_start = new AssignmentExpression(_start.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0));
            AssignmentExpression change_end = new AssignmentExpression(_end.clone(), AssignmentOperator.NORMAL, OneArith.minusOne(temp_length_id.clone()));
            if_c_st.addStatement(new ExpressionStatement(change_start));
            if_c_st.addStatement(new ExpressionStatement(change_end));

            List<Specifier> null_pointer_specs = new ArrayList<Specifier>(Arrays.asList(Specifier.VOID, PointerSpecifier.UNQUALIFIED));
            List<Specifier> pointer_specs = new ArrayList<Specifier>(Arrays.asList(spec, PointerSpecifier.UNQUALIFIED));
            // c_st.addDeclaration(createDecls.createVariableDeclaration(pointer_specs, temp, arr.clone()));
            if_c_st.addDeclaration(createDecls.createVariableDeclaration(pointer_specs, temp, new Typecast(null_pointer_specs, new IntegerLiteral(0))));


            IDExpression iter_id = new NameID("ss_i");
            Expression condition = new BinaryExpression(iter_id.clone(), BinaryOperator.COMPARE_LT, for_end.clone());
            Expression i_stride = new IntegerLiteral(1);

            Expression step = new AssignmentExpression(iter_id.clone(), AssignmentOperator.ADD, i_stride);
            Declaration init_decl = createDecls.createVariableDeclaration(Specifier.INT, iter_id.clone(), for_start.clone());
            Statement init_stmt = new DeclarationStatement(init_decl);
            // Statement init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, for_start.clone()));
            
            CompoundStatement for_body = new CompoundStatement();
            for_loop = new ForLoop(init_stmt, condition, step, for_body);
            
            if(pattern_ind[0]!=-1) {
                for_start = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(0));
                for_end = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(1));

                i_stride.swapWith(OneArith.addOne(new BinaryExpression(for_end, BinaryOperator.SUBTRACT, for_start)));

                iter_id = new NameID("ss_j");
                condition = new BinaryExpression(iter_id.clone(), BinaryOperator.COMPARE_LT, for_end.clone());
                step = new UnaryExpression(UnaryOperator.POST_INCREMENT, iter_id.clone());

                init_decl = createDecls.createVariableDeclaration(Specifier.INT, iter_id.clone(), for_start.clone());
                init_stmt = new DeclarationStatement(init_decl);
                // init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, for_start.clone()));
                CompoundStatement new_for_body = new CompoundStatement();
                ForLoop new_for_loop = new ForLoop(init_stmt, condition, step, new_for_body);
                for_body.addStatement(new_for_loop);
                for_body = new_for_body;
            }

            Pair<Expression, Expression> expr_nstride = createStridedExpr(true);
            new_stride = expr_nstride.getSecond();
            // AssignmentExpression change_stride = new AssignmentExpression(_stride, AssignmentOperator.NORMAL, new_stride.clone());


            
            
            Expression _left = new ArrayAccess(temp.clone(), expr_nstride.getFirst());
            Expression _right = new ArrayAccess(arr.clone(), createStridedExpr(false).getFirst());

            Expression body_expr = null;
            if((copy_direction==DpuCopy.DPU_COPY_TO) || (copy_direction==DpuCopy.DPU_BROADCAST_TO))
                body_expr = new AssignmentExpression(_left, AssignmentOperator.NORMAL, _right);
            else
                body_expr = new AssignmentExpression(_right, AssignmentOperator.NORMAL, _left);

            for_body.addStatement(new ExpressionStatement(body_expr));

            FunctionCall modify_dpu_copy_call = dpu_copy_call.clone();
            IRTools.replaceAll(modify_dpu_copy_call, arr, temp);

            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                
                if((new_stride!=null)&&(!new_stride.equals(new IntegerLiteral(1))))  {
                    IDExpression new_stride_name = new NameID(gpu_arr.getName()+"_stride_"+Integer.toString(arg_ind));
                    VariableDeclaration stride_decl = createDecls.createVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG)), new_stride_name.clone(), new_stride);
                    if_c_st.addDeclaration(stride_decl);

                    FunctionCall nfunc_call = HostOneTBtoDPU.createDpuCopy(new_stride_name, new IntegerLiteral(0), new_stride_name, null, null, ArgType.NORMAL, null, copy_direction);
                    if_c_st.addStatement(new ExpressionStatement(nfunc_call));
                }

                if_c_st.addStatement(for_loop);
                if_c_st.addStatement(new ExpressionStatement(modify_dpu_copy_call));
            }
            else {
                if_c_st.addStatement(new ExpressionStatement(modify_dpu_copy_call));
                if_c_st.addStatement(for_loop);
            }
            
            if_c_st.addStatement(new ExpressionStatement(createDecls.createMemFree(temp.clone()) ));
            IfStatement if_stmt = new IfStatement(cond, if_c_st, new ExpressionStatement(dpu_copy_call.clone()) );
            c_st.addStatement(if_stmt);
            // return for_loop;
        }
        
        return new Pair(c_st, new_stride);

    }

    static Pair<Statement, Expression> generateGeneralisedColAccCopy(Specifier spec, Expression arr, IDExpression gpu_arr, IDExpression temp, DpuCopy copy_direction, FunctionCall dpu_copy_call, Expression fstride, Literal copy_start) {
        CompoundStatement c_st = new CompoundStatement();
        Expression new_stride = null;
        IDExpression copy_id = new NameID("copy_i");
        IDExpression se_id = new NameID("ss_SE");

        String pre = "";
        if(copy_direction==DpuCopy.DPU_COPY_FROM || copy_direction==DpuCopy.DPU_BROADCAST_FROM) {
            pre = "p_";    
        }

        Expression _end = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_end"), copy_id.clone());
        Expression _start = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_offset"), copy_id.clone());
        Expression _start_cpu = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_cpu_offset"), copy_id.clone());
        Expression _stride = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_stride"), copy_id.clone());
        Expression _fstride = new NameID("f_"+pre+gpu_arr.getName()+"_stride");

        Expression _length = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_length"), copy_id.clone());
        Expression _flag = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_flag"), copy_id.clone());

            // System.out.println("Arg Type: "+spec);
            ForLoop for_loop = null;
            Expression dpu_copy_stride = new IntegerLiteral(1);
            // ArrayAccess for_start = new ArrayAccess(s_id[pattern_ind[1]].clone(), new IntegerLiteral(0));
            // ArrayAccess for_end = new ArrayAccess(s_id[pattern_ind[1]].clone(), new IntegerLiteral(1));

            IDExpression temp_length_id = new NameID(temp.getName()+"_length");
            IDExpression temp_ind = new NameID("temp_ind");
            // VariableDeclaration v_decl = createDecls.createVariableDeclaration(Specifier.INT, temp_length_id, tempBuffLength());
            
            // Expression _end = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_end"), new IntegerLiteral(arg_ind));
            // Expression _start = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_offset"), new IntegerLiteral(arg_ind));
            // Expression _stride = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_stride"), new IntegerLiteral(arg_ind));

            Expression length_row = OneArith.addOne(new BinaryExpression(_end.clone(), BinaryOperator.SUBTRACT, _start.clone()));
            Expression cond = null;
            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                cond = new BinaryExpression(length_row, BinaryOperator.COMPARE_GT, temp_length_id.clone());
            }
            else {
                // cond = new BinaryExpression(_stride.clone(), BinaryOperator.COMPARE_NE, new IntegerLiteral(1));
                cond = new BinaryExpression(_flag.clone(), BinaryOperator.BITWISE_AND, new IntegerLiteral(4));
            }

            CompoundStatement if_c_st = new CompoundStatement();
            CompoundStatement start_if_c_st = if_c_st;

            new_stride = fstride;

            

            Expression for_start = null;
            Expression for_end = null;


            Expression modify_length_expr = OneArith.addOne(new BinaryExpression(createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(3)), BinaryOperator.SUBTRACT, createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(1)) ));
            
            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                modify_length_expr = new BinaryExpression(modify_length_expr, BinaryOperator.MULTIPLY, _stride.clone());

            }
            else {
                modify_length_expr = _length;
            }

            Expression modify_length = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, modify_length_expr);
            Expression modify_length_align = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, KernelRepr.make8ByteAligned(temp_length_id.clone(), createDecls.createSizeofExpression(new NameID(spec.toString()))));
            // Expression modify_length_align = new AssignmentExpression(temp_length_id.clone(), AssignmentOperator.NORMAL, KernelRepr.make8ByteAligned(temp_length_id.clone(), element_size_id.clone()));
            
            VariableDeclaration t_ind_decl = createDecls.createVariableDeclaration(Specifier.INT, temp_ind, new IntegerLiteral(0));
            if_c_st.addDeclaration(t_ind_decl);

            AssignmentExpression change_start = new AssignmentExpression(_start.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0));
            AssignmentExpression change_end = new AssignmentExpression(_end.clone(), AssignmentOperator.NORMAL, OneArith.minusOne(temp_length_id.clone()));
            
            // if_c_st.addStatement(new ExpressionStatement(modify_length)); 

            if(copy_direction==DpuCopy.DPU_COPY_TO|| copy_direction==DpuCopy.DPU_BROADCAST_TO) {

                if_c_st.addStatement(new ExpressionStatement(change_start));
                if_c_st.addStatement(new ExpressionStatement(change_end));
                if_c_st.addStatement(new ExpressionStatement(new AssignmentExpression(_flag, AssignmentOperator.BITWISE_INCLUSIVE_OR, new IntegerLiteral(4)) ));

            }

            if(arr == null) {
                CompoundStatement else_comp_stmt = new CompoundStatement();
                if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                // else_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(_stride.clone(), AssignmentOperator.NORMAL, _fstride.clone()) ));
                    else_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(_stride.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(1)) ));
                }
                
                

                IfStatement if_stmt = new IfStatement(cond, start_if_c_st, else_comp_stmt);
                // c_st.addStatement(if_stmt);
                // return for_loop;
            
                return new Pair(if_stmt, new_stride);
            }

            VariableDeclaration v_decl = createDecls.createVariableDeclaration(Specifier.INT, temp_length_id, HostOneTBtoDPU.createSizeofRectange());
            c_st.addDeclaration(v_decl);
            
            if_c_st.addStatement(new ExpressionStatement(modify_length_align));
            


            List<Specifier> null_pointer_specs = new ArrayList<Specifier>(Arrays.asList(Specifier.VOID, PointerSpecifier.UNQUALIFIED));
            List<Specifier> pointer_specs = new ArrayList<Specifier>(Arrays.asList(spec, PointerSpecifier.UNQUALIFIED));
            // c_st.addDeclaration(createDecls.createVariableDeclaration(pointer_specs, temp, arr.clone()));
            Declaration null_pointer_decl = createDecls.createVariableDeclaration(pointer_specs, temp, new Typecast(null_pointer_specs, new IntegerLiteral(0)));

            if_c_st.addDeclaration(null_pointer_decl);
            // if(copy_start!=null) {
            //     Expression check_copy_ind_cond = new BinaryExpression(copy_id.clone(), BinaryOperator.COMPARE_GT, copy_start.clone());
            //     IfStatement check_copy_ind = new IfStatement(check_copy_ind_cond, new CompoundStatement());
            //     if_c_st.addStatement(check_copy_ind);
            //     if_c_st = (CompoundStatement)check_copy_ind.getThenStatement();
            // }

            ArrayAccess flag_id_acc = new ArrayAccess(new NameID(pre+gpu_arr.getName()+"_flag"), copy_id.clone());
            int flag_check_val = 0;
            if(copy_direction == DpuCopy.DPU_COPY_TO) {
                flag_check_val = 1;
            }
            else {
                flag_check_val = 2;
            }

            Expression check_rw_cond = new BinaryExpression(flag_id_acc.clone(), BinaryOperator.BITWISE_AND, new IntegerLiteral(flag_check_val));
            IfStatement check_rw = new IfStatement(check_rw_cond, new CompoundStatement());
            if_c_st.addStatement(check_rw);
            if_c_st = (CompoundStatement)check_rw.getThenStatement();

            if(copy_direction == DpuCopy.DPU_COPY_TO || copy_direction == DpuCopy.DPU_BROADCAST_TO) {
                for_start = createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(1));
                for_end = createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(3));
            }

            else {
                VariableDeclaration ij_decl = createDecls.createVariableDeclaration(Specifier.INT, new NameID("ss_i_start"), new BinaryExpression(_start_cpu, BinaryOperator.DIVIDE, fstride.clone()));
                if_c_st.addDeclaration(ij_decl);

                Expression i_limit = new BinaryExpression(_length.clone(), BinaryOperator.DIVIDE, _stride.clone());
                ij_decl = createDecls.createVariableDeclaration(Specifier.INT, new NameID("ss_i_limit"), new BinaryExpression(new NameID("ss_i_start"), BinaryOperator.ADD, i_limit));
                if_c_st.addDeclaration(ij_decl);
                
                ij_decl = createDecls.createVariableDeclaration(Specifier.INT, new NameID("ss_j_start"), new BinaryExpression(_start_cpu.clone(), BinaryOperator.MODULUS, fstride.clone()));
                if_c_st.addDeclaration(ij_decl);

                ij_decl = createDecls.createVariableDeclaration(Specifier.INT, new NameID("ss_j_limit"), new BinaryExpression(_stride.clone(), BinaryOperator.ADD, new NameID("ss_j_start")));
                if_c_st.addDeclaration(ij_decl);

                String code_print_ij = "//printf(\"ss_i:%d-%d, ss_j:%d-%d\\n\",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);";
                if_c_st.addStatement(createDecls.createCodeAnnotStmt(code_print_ij));   

                for_start = new NameID("ss_i_start");
                for_end = new NameID("ss_i_limit");
            }

            
            if_c_st.addStatement(new ExpressionStatement(createMallocTemp(spec, temp, temp_length_id.clone())) );

            IDExpression iter_id = new NameID("ss_i");
            BinaryOperator iter_op = BinaryOperator.COMPARE_LT;
            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                iter_op = BinaryOperator.COMPARE_LE;
            }
            Expression condition = new BinaryExpression(iter_id.clone(), iter_op, for_end);
            Expression i_stride = new IntegerLiteral(1);

            Expression step = new AssignmentExpression(iter_id.clone(), AssignmentOperator.ADD, i_stride);
            Declaration init_decl = createDecls.createVariableDeclarationNC(Specifier.INT, iter_id.clone(), for_start);
            Statement init_stmt = new DeclarationStatement(init_decl);
            // Statement init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, for_start.clone()));
            
            CompoundStatement for_body = new CompoundStatement();
            for_loop = new ForLoop(init_stmt, condition, step, for_body);
            
            // if(pattern_ind[0]!=-1) {
                // for_start = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(0));
                // for_end = new ArrayAccess(s_id[pattern_ind[0]].clone(), new IntegerLiteral(1));
                if(copy_direction == DpuCopy.DPU_COPY_TO || copy_direction == DpuCopy.DPU_BROADCAST_TO) {
                    for_start = createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(0));
                    for_end = createDecls.createArrayAccess(se_id.clone(), copy_id.clone(), new IntegerLiteral(2));
                }
                else {
                    for_start = new NameID("ss_j_start");
                    for_end = new NameID("ss_j_limit");
                }

                dpu_copy_stride = OneArith.addOne(new BinaryExpression(for_end, BinaryOperator.SUBTRACT, for_start));

                // i_stride.swapWith(OneArith.addOne(new BinaryExpression(for_end, BinaryOperator.SUBTRACT, for_start)));
                // i_stride.swapWith(_stride);

                iter_id = new NameID("ss_j");
                condition = new BinaryExpression(iter_id.clone(), iter_op, for_end.clone());
                step = new UnaryExpression(UnaryOperator.POST_INCREMENT, iter_id.clone());

                init_decl = createDecls.createVariableDeclaration(Specifier.INT, iter_id.clone(), for_start.clone());
                init_stmt = new DeclarationStatement(init_decl);
                // init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, for_start.clone()));
                CompoundStatement new_for_body = new CompoundStatement();
                ForLoop new_for_loop = new ForLoop(init_stmt, condition, step, new_for_body);
                for_body.addStatement(new_for_loop);
                for_body = new_for_body;
            // }

            
            // AssignmentExpression change_stride = new AssignmentExpression(_stride, AssignmentOperator.NORMAL, new_stride.clone());


            
            
            // Expression _left = new ArrayAccess(temp.clone(), expr_nstride.getFirst());
            CompoundStatement else_comp_stmt = new CompoundStatement();

            if(arr!=null) {


                Expression temp_ind_inc = new UnaryExpression(UnaryOperator.POST_INCREMENT, temp_ind.clone());
                Expression _left = new ArrayAccess(temp.clone(), temp_ind_inc);

                BinaryExpression right_ind = new BinaryExpression(fstride.clone(), BinaryOperator.MULTIPLY, new NameID("ss_i"));
                right_ind = new BinaryExpression(right_ind, BinaryOperator.ADD, new NameID("ss_j"));

                Expression _right = new ArrayAccess(arr.clone(), right_ind);

                Expression body_expr = null;
                if((copy_direction==DpuCopy.DPU_COPY_TO) || (copy_direction==DpuCopy.DPU_BROADCAST_TO)) {
                    body_expr = new AssignmentExpression(_left, AssignmentOperator.NORMAL, _right);
                    for_body.addStatement(new ExpressionStatement(body_expr));
                }
                else {
                    body_expr = new AssignmentExpression(_right, AssignmentOperator.NORMAL, _left);
                    for_body.addStatement(new ExpressionStatement(body_expr));
                    Expression overflow_check_cond = new BinaryExpression(temp_ind.clone(), BinaryOperator.COMPARE_GE, _length.clone());
                    IfStatement overflow_check = new IfStatement(overflow_check_cond, new BreakStatement());
                    for_body.addStatement(overflow_check);
                }

                FunctionCall modify_dpu_copy_call = dpu_copy_call.clone();
                modify_dpu_copy_call.setArgument(3, temp.clone());
                IRTools.replaceAll(modify_dpu_copy_call, arr, temp);


                if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                
                // if((new_stride!=null)&&(!new_stride.equals(new IntegerLiteral(1))))  {
                //     IDExpression new_stride_name = new NameID(gpu_arr.getName()+"_stride");
                //     // VariableDeclaration stride_decl = createDecls.createVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG)), new_stride_name.clone(), new_stride);
                //     VariableDeclaration stride_decl = createDecls.createVariableDeclaration(new ArrayList<Specifier>(Arrays.asList(Specifier.LONG)), new_stride_name.clone(), dpu_copy_stride);
                //     if_c_st.addDeclaration(stride_decl);

                //     // if(!dpu_copy_stride.equals(new IntegerLiteral(1))) {
                //         // FunctionCall nfunc_call = HostOneTBtoDPU.createDpuCopy(new_stride_name, new IntegerLiteral(0), new_stride_name, null, ArgType.NORMAL, null, copy_direction);
                //         // if_c_st.addStatement(new ExpressionStatement(nfunc_call));
                //     // }
                // }

                    if_c_st.addStatement(for_loop);
                    if_c_st.addStatement(new ExpressionStatement(modify_dpu_copy_call));
                }
                else {
                    if_c_st.addStatement(new ExpressionStatement(modify_dpu_copy_call));
                    if_c_st.addStatement(for_loop);
                }

                if_c_st.addStatement(new ExpressionStatement(createDecls.createMemFree(temp.clone()) ));
            

                FunctionCall else_dpu_call = dpu_copy_call.clone();
                IDExpression copy_size = null;
                if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO)
                    copy_size = new NameID("copy_to_size");
                else
                    copy_size = new NameID("copy_from_size");

                Expression else_copy_size = else_dpu_call.getArgument(4);
                else_dpu_call.setArgument(4, copy_size);


                else_comp_stmt.addDeclaration(createDecls.createVariableDeclaration(Specifier.INT, copy_size.clone(), else_copy_size));
            
                CompoundStatement else_copy_stmt = new CompoundStatement();

                if(copy_direction==DpuCopy.DPU_COPY_FROM || copy_direction==DpuCopy.DPU_BROADCAST_FROM) {

                    else_copy_stmt.addDeclaration(null_pointer_decl.clone());
                    else_copy_stmt.addDeclaration(v_decl.clone());

                    Expression overflow_check_cond = new BinaryExpression(copy_size.clone(), BinaryOperator.COMPARE_GT, _length.clone());
                    Expression malloc_temp = createDecls.createMallocVariableAssign(spec, temp, new BinaryExpression(copy_size.clone(), BinaryOperator.SUBTRACT, _length.clone()));
                    CompoundStatement overflow_check_c_st = new CompoundStatement();
                    overflow_check_c_st.addStatement(new ExpressionStatement(malloc_temp));

                    CompoundStatement overflow_check_loop_body = new CompoundStatement();
                    IDExpression loop_id = new NameID("ss_i");
                    Statement loop_init = new DeclarationStatement(createDecls.createVariableDeclaration(Specifier.INT, loop_id.clone(), new BinaryExpression(_start_cpu.clone(), BinaryOperator.ADD, _length.clone())));
                    Expression loop_cond = new BinaryExpression(loop_id.clone(), BinaryOperator.COMPARE_LT, copy_size.clone());
                    Expression loop_inc = new UnaryExpression(UnaryOperator.POST_INCREMENT, loop_id.clone());

                    Expression cp_left = new ArrayAccess(arr.clone(), new NameID("ss_i"));
                    Expression cp_right = new ArrayAccess(temp.clone(), new UnaryExpression(UnaryOperator.POST_INCREMENT, temp_ind.clone()));
                    overflow_check_loop_body.addStatement(new ExpressionStatement(new AssignmentExpression(cp_right, AssignmentOperator.NORMAL, cp_left)));

                    ForLoop overflow_check_loop = new ForLoop(loop_init, loop_cond, loop_inc, overflow_check_loop_body);
                    overflow_check_c_st.addStatement(overflow_check_loop);

                    IfStatement overflow_check = new IfStatement(overflow_check_cond, overflow_check_c_st);
                    else_copy_stmt.addStatement(overflow_check);
                }
                
                else_copy_stmt.addStatement(new ExpressionStatement(else_dpu_call));
                else_copy_stmt.addDeclaration(t_ind_decl.clone());

                if(copy_direction==DpuCopy.DPU_COPY_FROM || copy_direction==DpuCopy.DPU_BROADCAST_FROM) {
                    Expression overflow_check_cond = new BinaryExpression(copy_size.clone(), BinaryOperator.COMPARE_GT, _length.clone());
                    CompoundStatement overflow_check_c_st = new CompoundStatement();
                    overflow_check_c_st.addStatement(new ExpressionStatement(new AssignmentExpression(temp_ind.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0)) ));

                    CompoundStatement overflow_check_loop_body = new CompoundStatement();
                    IDExpression loop_id = new NameID("ss_i");
                    Statement loop_init = new DeclarationStatement(createDecls.createVariableDeclaration(Specifier.INT, loop_id.clone(), new BinaryExpression(_start_cpu.clone(), BinaryOperator.ADD, _length.clone())));
                    Expression loop_cond = new BinaryExpression(loop_id.clone(), BinaryOperator.COMPARE_LT, copy_size.clone());
                    Expression loop_inc = new UnaryExpression(UnaryOperator.POST_INCREMENT, loop_id.clone());
                    Expression cp_left = new ArrayAccess(arr.clone(), loop_id.clone());
                    Expression cp_right = new ArrayAccess(temp.clone(), new UnaryExpression(UnaryOperator.POST_INCREMENT, temp_ind.clone()));
                    
                    overflow_check_loop_body.addStatement(new ExpressionStatement(new AssignmentExpression(cp_left, AssignmentOperator.NORMAL, cp_right)));

                    ForLoop overflow_check_loop = new ForLoop(loop_init, loop_cond, loop_inc, overflow_check_loop_body);
                    overflow_check_c_st.addStatement(overflow_check_loop);

                    overflow_check_c_st.addStatement(new ExpressionStatement(createDecls.createMemFree(temp.clone()) ));

                    IfStatement overflow_check = new IfStatement(overflow_check_cond, overflow_check_c_st);
                    else_copy_stmt.addStatement(overflow_check);
                }

                check_rw = new IfStatement(check_rw_cond.clone(), else_copy_stmt);
                else_comp_stmt.addStatement(check_rw);
            }
            
            // if(copy_direction==DpuCopy.DPU_COPY_TO|| copy_direction==DpuCopy.DPU_BROADCAST_TO)
            //     if_c_st.addStatement(new ExpressionStatement(new AssignmentExpression(_flag, AssignmentOperator.BITWISE_INCLUSIVE_OR, new IntegerLiteral(4)) ));
            

            

            // if(copy_start!=null) {
            //     Expression check_copy_ind_cond = new BinaryExpression(copy_id.clone(), BinaryOperator.COMPARE_GT, copy_start.clone());
            //     IfStatement check_copy_ind = new IfStatement(check_copy_ind_cond, else_comp_stmt);
            //     else_comp_stmt = check_copy_ind;
            // }

            
            
            
            if(copy_direction==DpuCopy.DPU_COPY_TO || copy_direction==DpuCopy.DPU_BROADCAST_TO) {
                // else_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(_stride.clone(), AssignmentOperator.NORMAL, _fstride.clone()) ));
                else_comp_stmt.addStatement(new ExpressionStatement(new AssignmentExpression(_stride.clone(), AssignmentOperator.NORMAL, new IntegerLiteral(1)) ));
            }
            
            

            IfStatement if_stmt = new IfStatement(cond, start_if_c_st, else_comp_stmt);
            // c_st.addStatement(if_stmt);
            // return for_loop;
        
        return new Pair(if_stmt, new_stride);

    }

    void findExprStrides(Expression expr, Primitive primitive, KernelRepr _kernel) {
        // Strides ind_coefficients = new Strides();
        // Expression analyse_expr = Symbolic.simplify(expr.clone());
        Expression analyse_expr = expr.clone();
        IDExpression div_group_grid = _kernel.get_partition_dim_grid();

        System.out.println("Analysed Expression = "+analyse_expr);
        Triple<Integer, Integer, Boolean> res_coefficients = AnalyseArrayExpression.analyseIndex(analyse_expr, this, 0, primitive, _kernel);
        if(res_coefficients.getFirst()==-1) {
            this.is_strided = -1;
            System.out.println("[ERROR] check array index for mul variables" + analyse_expr);
        }
        else {
            this.set_offset(analyse_expr);
            this.simplify_all();
            this.get_overall_stride(div_group_grid);
            this.print();
        }    
    }
    



}
class TranslationVariableNames {

    static IDExpression getIteratorVariableHost(long dim) {
        IDExpression iter_id = null;
        if(dim>9) {
            System.out.println("ERROR: dimension provided is a value not supported");
        }

//         String[] id_names = new String[]{"ss_tx", "ss_ty", "ss_tz", "ss_gx", "ss_gy", "ss_gz", "ss_bx", "ss_by", "ss_bz"};
        String[] id_names = new String[]{"tx", "ty", "tz", "gx", "gy", "gz", "bx", "by", "bz"};
        iter_id = new NameID(id_names[(int)dim]);

           
        return iter_id;
    }
    static IDExpression getIteratorVariableHost(int dim, String pre, String post) {
        IDExpression iter_id = null;
        if(dim>8) {
            System.out.println("ERROR: dimension provided is a value not supported");
        }
        
//             String[] id_names = new String[]{"ss_tx", "ss_ty", "ss_tz", "ss_gx", "ss_gy", "ss_gz", "ss_bx", "ss_by", "ss_bz"};
        String[] id_names = new String[]{"tx", "ty", "tz", "gx", "gy", "gz", "bx", "by", "bz"};
        
        if(post==null)
            iter_id = new NameID(pre+id_names[dim]);
        else
            iter_id = new NameID(pre+id_names[dim]+post);
        
        return iter_id;
    }
    
    static Expression getIteratorLimit(long dim) { // gives the start of tx,ty,tz for a block
//         String[] limit_names = new String[]{"ss_tx_end", "ss_ty_end", "ss_tz_end", "ss_gx_end", "ss_gy_end", "ss_gz_end"};

        if(dim>3)
            return null;
        
        ArrayAccess limit_name = new ArrayAccess(new NameID("t_td_se"), new IntegerLiteral(dim));
        Expression limit_limit = new BinaryExpression(limit_name, AccessOperator.MEMBER_ACCESS, new NameID("end"));
        
        return limit_limit;
        
    }
    
    static Expression getIteratorStart(long dim) {     // gives the start of tx,ty,tz for a block
//         String[] limit_start_names = new String[]{"ss_tx_start", "ss_ty_start", "ss_tz_start", "ss_gx_start", "ss_gy_start", "ss_gz_start"};
        if(dim>3)
            return null;
        
        ArrayAccess limit_name = new ArrayAccess(new NameID("t_td_se"), new IntegerLiteral(dim));
        Expression limit_start = new BinaryExpression(limit_name, AccessOperator.MEMBER_ACCESS, new NameID("start"));
        
        return limit_start;
        
    }
}

class WorkFuncRangeMap{
    // Map<IDExpression, Long> myMap = new HashMap<IDExpression, Long>() {{
    //     put(new NameID("txl"), null);
    //     put(new NameID("tyl"), null);
    //     put(new NameID("tzl"), null);

    //     put(new NameID("gxl"), null);
    //     put(new NameID("gyl"), null);
    //     put(new NameID("gzl"), null);

    //     put(new NameID("bxl"), null);
    //     put(new NameID("byl"), null);
    //     put(new NameID("bzl"), null);
    // }};

    Map<String, Pair<Expression,Expression>> myMap = new HashMap<String, Pair<Expression, Expression>>() {{
        put("txl", null);
        put("tyl", null);
        put("tzl", null);

        put("gxl", null);
        put("gyl", null);
        put("gzl", null);

        put("bxl", null);
        put("byl", null);
        put("bzl", null);
    }};

    Expression createArrayAcc(String arr_name, int ind) {
        ArrayAccess arr_acc = new ArrayAccess(new NameID(arr_name), new IntegerLiteral(ind));
        return arr_acc;
    }

    void fillNulls(String xdl) {
        Pair p = myMap.get(xdl);
        if(p==null) {
            myMap.put(xdl, new Pair(createArrayAcc(xdl, 0), createArrayAcc(xdl, 1)));
        }
    }

    Expression createmMFuncCall(int is_min, Expression expr1, Expression expr2) {
        if(is_min==2) 
            return expr2.clone();

        IDExpression[] func_names = {new NameID("ss_max"), new NameID("ss_min")};
        return new FunctionCall(func_names[is_min], expr1.clone(), expr2.clone());


    }

    boolean addMapKeys(List<String> convert_ones) {
        for(String str: convert_ones) {
            if(!myMap.containsKey(str)){
                myMap.put(str, null);
            }
        }
        return true;
    }

    boolean setMap(String xdl, int ind, Expression expr, int is_min) {
        if(myMap.containsKey(xdl)==false)
            return false;
        fillNulls(xdl);
        
        
        if(ind==0) {
            Expression cur_expr = myMap.get(xdl).getFirst();
            expr = createmMFuncCall(is_min, cur_expr, expr);
            myMap.get(xdl).setFirst(expr);
        }
        else if(ind==1){
            Expression cur_expr = myMap.get(xdl).getSecond();
            expr = createmMFuncCall(is_min, cur_expr, expr);
            myMap.get(xdl).setSecond(expr);
        }
        return true;
    }

    static Expression modifyToKexpr(KernelRepr _kernel, Primitive k_primitive, Expression expr) {
        // if(expr instanceof IDExpression) {
        //     IDExpression id_expr = (IDExpression) expr;
        //     if(k_primitive.total_parameter_ind.containsKey(id_expr)) {
        //         int index = k_primitive.total_parameter_ind.get(id_expr);
        //         return _kernel.argument_list.get(index).getSecond();
        //     }
        //     return expr.clone();
        // }
        Expression dup_expr = expr.clone();
        // DFIterator<IDExpression> df_iterator = new DFIterator<IDExpression>(dup_expr, IDExpression.class);

        // while(df_iterator.hasNext()) {
        //     IDExpression id_expr = df_iterator.next();
        //     if(k_primitive.total_parameter_ind.containsKey(id_expr)) {
        //         int index = k_primitive.total_parameter_ind.get(id_expr);
        //         id_expr.swapWith(_kernel.argument_list.get(index).getSecond());
        //     }
        // }
        return dup_expr;
    }

    static Statement modifyToKexpr(KernelRepr _kernel, Primitive k_primitive, Statement expr) {
        Statement dup_expr = expr.clone();
        DFIterator<IDExpression> df_iterator = new DFIterator<IDExpression>(dup_expr, IDExpression.class);

        while(df_iterator.hasNext()) {
            IDExpression id_expr = df_iterator.next();
            if(k_primitive.total_parameter_ind.containsKey(id_expr)) {
                int index = k_primitive.total_parameter_ind.get(id_expr);
                id_expr.swapWith(_kernel.argument_list.get(index).getSecond());
            }
        }
        return dup_expr;
    }

    List<Declaration> createWorkFuncDecls(KernelRepr _kernel, Primitive k_primitive) {
        List<Declaration> wf_decls = new ArrayList<Declaration>();

        for(Map.Entry<String, Pair<Expression, Expression>> entry : myMap.entrySet()) {
            if(entry.getValue()==null)
                continue;
            // Expression l_limit = modifyToKexpr(_kernel, k_primitive, entry.getValue().getFirst());
            // Expression r_limit = modifyToKexpr(_kernel, k_primitive, entry.getValue().getSecond());

            Expression l_limit = entry.getValue().getFirst();
            Expression r_limit = entry.getValue().getSecond();

            Declaration decl = createDecls.createArrayVariableDeclaration(Specifier.LONG, new NameID("_"+entry.getKey()), l_limit, r_limit).getDeclaration();
            wf_decls.add(decl);
        }

        return wf_decls;
    }

    void modify_expr(Traversable trav){
        List<String> convert_ones = new ArrayList(Arrays.asList("txl", "tyl", "tzl", "bxl", "byl", "bzl", "gxl", "gyl", "gzl"));

        DFIterator<Statement> df_st_iterator = new DFIterator<Statement>(trav, Statement.class);
        while(df_st_iterator.hasNext()) {
            Statement st = df_st_iterator.next();

            if(st instanceof ExpressionStatement) {
                DFIterator<ArrayAccess> df_iterator = new DFIterator<ArrayAccess>(st, ArrayAccess.class);

                while(df_iterator.hasNext()) {
                    ArrayAccess arr_acc = df_iterator.next();

                    if(arr_acc.getParent() instanceof Initializer)
                        continue;
                    if(!(arr_acc.getArrayName() instanceof IDExpression))
                        continue;

                    IDExpression _id = (IDExpression)arr_acc.getArrayName();

                    if(convert_ones.indexOf(_id.getName())!=-1) {
                        if(myMap.get(_id.getName())!=null) {
                            _id.swapWith(new NameID("_"+_id.getName()));
                        }
                    }
                }
            }
        }
        

    }

    void print() {
        for(Map.Entry<String, Pair<Expression, Expression>> entry : myMap.entrySet()) {
            if(entry.getValue()==null)
                continue;
            System.out.print(entry.getKey()+" ");
            entry.getValue().print();
        }
    }

}

class createDecls{
    static Declaration createProcedure(Object o_spec, String s_proc_id, List<String> s_param_int, List<String> s_param_pointer_int, List<String> s_param_pointer_Point, Statement body) {
    
        List<Specifier> p_spec = null;
        if(o_spec instanceof Specifier) {
            p_spec = new ArrayList<Specifier>(Arrays.asList((Specifier)o_spec));
        }
        else if(o_spec instanceof List) {
            p_spec = (ArrayList<Specifier>)o_spec;
        }
        IDExpression proc_id = new NameID(s_proc_id);
        
        
        List<VariableDeclaration> params = new ArrayList<VariableDeclaration>();
        
        toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), s_param_int, null);
        toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(Specifier.INT, PointerSpecifier.UNQUALIFIED)), s_param_pointer_int, null);
        toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(new UserSpecifier(new NameID("Point")), PointerSpecifier.UNQUALIFIED)), s_param_pointer_Point, null);
        
        ProcedureDeclarator p_decl = new ProcedureDeclarator( proc_id, params);
        
        CompoundStatement c_body = null;
        if(body instanceof CompoundStatement)
            c_body = (CompoundStatement)body;
        else {
            c_body = new CompoundStatement();
            c_body.addStatement(body);
        }
        
        return new Procedure(p_spec, p_decl, c_body);
    }

    static Declaration createProcedure(List<Specifier> p_spec, List<Specifier> trail_spec, String s_proc_id, List<String> s_param_int_trail, List<String> s_param_pointer_int, List<String> s_param_int, Statement body) {
    
        IDExpression proc_id = new NameID(s_proc_id);
        
        
        List<VariableDeclaration> params = new ArrayList<VariableDeclaration>();
        
        if(trail_spec==null)
            toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), s_param_int_trail, null);
        else 
            toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), s_param_int_trail, trail_spec);

        toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(Specifier.LONG, Specifier.LONG, PointerSpecifier.UNQUALIFIED)), s_param_pointer_int, null);
        toVariableDeclaration(params, new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), s_param_int, null);        
        ProcedureDeclarator p_decl = new ProcedureDeclarator(proc_id, params);
        
        CompoundStatement c_body = null;
        if(body instanceof CompoundStatement)
            c_body = (CompoundStatement)body;
        else {
            c_body = new CompoundStatement();
            c_body.addStatement(body);
        }
        
        return new Procedure(p_spec, p_decl, c_body);
    }
    
    List<IDExpression> toIDExpression(List<String> string_vars) {
        List<IDExpression> id_vars = new ArrayList<IDExpression>();
        if(string_vars!=null) {
            for(String var: string_vars) {
                id_vars.add(new NameID(var));
            }
        }
        return id_vars;
    }
    
    static void toVariableDeclaration(List<VariableDeclaration> decl_vars, List<Specifier> spec, List<String> string_vars, List<Specifier> trail_specs) {
        if(string_vars!=null) {
            for(String var: string_vars) {
                decl_vars.add(createVariableDeclaration(spec, new NameID(var), null, trail_specs));
            }
        }
//         return decl_vars;
    }

    static VariableDeclaration createVariableDeclarationNC(Specifier spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator declr = new VariableDeclarator(iter_id);
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        if(init_expr != null) {
            Initializer init = new Initializer(init_expr);
            declr.setInitializer(init);
        }
        return v_decl;
    }

     static VariableDeclaration createVariableDeclarationNC(Specifier spec, IDExpression iter_id, Initializer init) {
        VariableDeclarator declr = new VariableDeclarator(iter_id);
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        if(init != null) {
            declr.setInitializer(init);
        }
        return v_decl;
    }
    
    static VariableDeclaration createVariableDeclaration(List<Specifier> spec, IDExpression iter_id, Expression init_expr, List<Specifier> trail_spec) {
        VariableDeclarator declr = null;
        if(trail_spec==null)
            declr = new VariableDeclarator(iter_id);
        else
            declr = new VariableDeclarator(iter_id, trail_spec);

        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        if(init_expr != null) {
            Initializer init = new Initializer(init_expr.clone());
            declr.setInitializer(init);
        }
        return v_decl;
    }

    static VariableDeclaration createVariableDeclaration(List<Specifier> spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator declr = new VariableDeclarator(iter_id);
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        if(init_expr != null) {
            Initializer init = new Initializer(init_expr.clone());
            declr.setInitializer(init);
        }
        return v_decl;
    }

    static VariableDeclaration createVariableDeclaration(Specifier spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator declr = new VariableDeclarator(iter_id);
        if(init_expr != null) {
            Initializer init = new Initializer(init_expr.clone());
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        return v_decl;
    }

        
    static DeclarationStatement createVariableDeclarationStmt(Specifier spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator declr = new VariableDeclarator(iter_id);
        if(init_expr != null) {
            Initializer init = new Initializer(init_expr.clone());
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        return init_iter;
    }

    static DeclarationStatement createVariableDeclarationStmt(List<Specifier> spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator declr = new VariableDeclarator(iter_id);
        if(init_expr != null) {
            Initializer init = new Initializer(init_expr.clone());
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        return init_iter;
    }

    static VariableDeclarator createArrayVariableDeclarator(Specifier spec, IDExpression iter_id, int size, Expression init_val) {
        ArraySpecifier arr_spec = ArraySpecifier.UNBOUNDED;
        if(size!=0)
            arr_spec = new ArraySpecifier(new IntegerLiteral(size));

        VariableDeclarator declr = new VariableDeclarator(iter_id, arr_spec);
        // VariableDeclaration v_decl = new VariableDeclaration(spec, declr);

        List<Expression> init_vals = new ArrayList<Expression>();
        if(init_val != null) {
            for(int i=0; i<size; i++) {
                init_vals.add(init_val.clone());
            }

            Initializer init = new Initializer(init_vals);
            declr.setInitializer(init);
        }

        return declr;
    }

    static Declaration createArrayVariableDeclaration(Specifier spec, Expression size, IDExpression iter_id) {
        ArraySpecifier arr_spec = ArraySpecifier.UNBOUNDED;
        if(size!=null)
            arr_spec = new ArraySpecifier(size);

        VariableDeclarator declr = new VariableDeclarator(iter_id, arr_spec);
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        return v_decl;
    }

    static Declaration createArrayVariableDeclaration(List<Specifier> specs, Expression size, IDExpression iter_id) {
        ArraySpecifier arr_spec = ArraySpecifier.UNBOUNDED;
        if(size!=null)
            arr_spec = new ArraySpecifier(size);

        VariableDeclarator declr = new VariableDeclarator(iter_id, arr_spec);
        VariableDeclaration v_decl = new VariableDeclaration(specs, declr);
        return v_decl;
    }
    
    static DeclarationStatement createArrayVariableDeclaration(Specifier spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator declr = new VariableDeclarator(iter_id, ArraySpecifier.UNBOUNDED);
        if(init_expr != null) {
            Initializer init = new Initializer(new ArrayList<Expression>() {{
                                    add(new IntegerLiteral(0));
                                    add(init_expr.clone());
                                    }}
                                    );
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        return init_iter;
    }
    
    static DeclarationStatement createArrayVariableDeclaration(Specifier spec, IDExpression iter_id, Expression init_expr1, Expression init_expr2) {
        VariableDeclarator declr = new VariableDeclarator(iter_id, ArraySpecifier.UNBOUNDED);
        if(init_expr1 != null && init_expr2!=null) {
            Initializer init = new Initializer(new ArrayList<Expression>() {{
                                    add(init_expr1.clone());
                                    add(init_expr2.clone());
                                    }}
                                    );
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        return init_iter;
    }
    
    static Declaration createArrayVariableDeclaration(Specifier spec, IDExpression iter_id, List<Expression> init_exprs) {
        VariableDeclarator declr = new VariableDeclarator(iter_id, ArraySpecifier.UNBOUNDED);
        if(init_exprs!=null) {
            
            
            Initializer init = new Initializer(init_exprs);
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        // DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        return v_decl;
    }

    static Declaration createArrayVariableDeclaration(List<Specifier> spec, IDExpression iter_id, List<Expression> init_exprs) {
        VariableDeclarator declr = new VariableDeclarator(iter_id, ArraySpecifier.UNBOUNDED);
        if(init_exprs!=null) {
            
            
            Initializer init = new Initializer(init_exprs);
            declr.setInitializer(init);
        }
            
//             Initializer init = new Initializer(getIteratorStart(dim));
        
        VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
        // DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        return v_decl;
    }

    static ArrayAccess createArrayAccess(IDExpression arr_id, Expression expr1, Expression expr2) {
        ArrayAccess arr_acc = new ArrayAccess(arr_id, new ArrayList<Expression>(Arrays.asList(expr1, expr2)));
        return arr_acc;
    }

    static Declaration createMallocVariableDeclaration(Specifier spec, IDExpression iter_id, Expression length) {
        VariableDeclarator temp_v_decl = new VariableDeclarator(PointerSpecifier.UNQUALIFIED, iter_id.clone());
        NestedDeclarator temp_nest_declr = new NestedDeclarator(temp_v_decl);
        
        Expression init_expr = new BinaryExpression(length.clone(), BinaryOperator.MULTIPLY, createSizeofExpression(new NameID(spec.toString())) );
        FunctionCall init_malloc = new FunctionCall(new NameID("malloc"), init_expr);
        temp_nest_declr.setInitializer(new Initializer(init_malloc));
        VariableDeclaration v_decl = new VariableDeclaration(spec, temp_nest_declr);
        return v_decl;
    }

    static Expression createMallocVariableAssign(Specifier spec, IDExpression iter_id, Expression length) {
        if(length==null)
            System.out.println("[ERROR] length=null");

        Expression init_expr = new BinaryExpression(length, BinaryOperator.MULTIPLY, createSizeofExpression(new NameID(spec.toString())) );
        FunctionCall init_malloc = new FunctionCall(new NameID("malloc"), init_expr);
        
        return new AssignmentExpression(iter_id.clone(), AssignmentOperator.NORMAL, init_malloc);
        
    }

    static Expression createMemFree(IDExpression _id) {
        return new FunctionCall(new NameID("free"), _id);
    }

    static Declaration createPointerDeclaration(Specifier spec, IDExpression iter_id, Expression init_expr) {
        VariableDeclarator temp_v_decl = new VariableDeclarator(PointerSpecifier.UNQUALIFIED, iter_id.clone());
        NestedDeclarator temp_nest_declr = new NestedDeclarator(temp_v_decl);
        
        if(init_expr!=null)
            temp_nest_declr.setInitializer(new Initializer(init_expr));

        VariableDeclaration v_decl = new VariableDeclaration(spec, temp_nest_declr);
        return v_decl;
    }

    static Declaration createProcedure(Specifier p_spec, String s_proc_id, boolean is_pointer, String param1, String param2, Statement body) {
    
        IDExpression proc_id = new NameID(s_proc_id);
    
    
        List<VariableDeclaration> params = new ArrayList<VariableDeclaration>();

        if(param1!="") {
            if(is_pointer) {
                List<Specifier> spec = new ArrayList<Specifier>(Arrays.asList(Specifier.INT, PointerSpecifier.UNQUALIFIED));
                VariableDeclarator declr = new VariableDeclarator(new NameID(param1));
                VariableDeclaration v_decl = new VariableDeclaration(spec, declr);
                params.add(v_decl);
            }
            else
                params.add(createVariableDeclaration(Specifier.INT, new NameID(param1), null));
        }
        if(param2!="") {
            params.add(createVariableDeclaration(Specifier.INT, new NameID(param2), null));
        }
    
        ProcedureDeclarator p_decl = new ProcedureDeclarator( proc_id, params);
    
    
        CompoundStatement c_body = null;

        if(body instanceof CompoundStatement)
            c_body = (CompoundStatement)body;
        else {
            c_body = new CompoundStatement();
            c_body.addStatement(body);
        }
    
        return new Procedure(new ArrayList<Specifier>(Arrays.asList(p_spec)), p_decl, c_body);
    }

    static FunctionCall createSizeofExpression(Expression expr) {
        return new FunctionCall(new NameID("sizeof"), expr);
    }
        
    
    static Statement createCodeAnnotStmt(String code) {
        CodeAnnotation c_annot = new CodeAnnotation(code);
        AnnotationStatement annot_stmt = new AnnotationStatement(c_annot);
        return annot_stmt;
    }
    
    static Declaration createCodeAnnotDecl(String code) {
        CodeAnnotation c_annot = new CodeAnnotation(code);
        AnnotationDeclaration annot_decl = new AnnotationDeclaration(c_annot);
        return annot_decl;
    }

}

class Primitive {
    HashSet<IDExpression> parameter;
    Map<IDExpression, Integer> parameter_map;
    Map<Integer, IDExpression> total_parameter_ind_map;
    Map<Integer, Specifier> total_parameter_spec_map;
    Map<IDExpression, Integer> total_parameter_ind;
//     HashSet<IDExpression> iterator_var;
    Map<IDExpression, RangeExpression> iterator_var;
    Map<IDExpression, Integer> iterator_var_map;
    List<IDExpression> global_buffers;
    List<IDExpression> local_buffers;
    List<IDExpression> scalars;
    
    int iterator_counter = 0;
//     HashSet<Expression> work_item_func;
    // String[] work_item_func_names = {"get_local_id", "get_global_id", "get_group_id"};
    // String[] work_item_variables = {"txl", "tyl", "tzl", "gxl", "gyl", "gzl", "bxl", "byl", "bzl"};
    List<String> work_item_func_names = new ArrayList<String>(Arrays.asList("get_local_id", "get_global_id", "get_group_id"));
    List<String> work_item_variables = new ArrayList<String>(Arrays.asList("txl", "tyl", "tzl", "gxl", "gyl", "gzl", "bxl", "byl", "bzl"));
    HashMap<IDExpression, Integer> work_item_func; 
    HashMap<Integer, IDExpression> work_item_map;
    
    Primitive() {
        parameter = new HashSet<IDExpression>();
        work_item_func = new HashMap<IDExpression, Integer>();
        work_item_map = new HashMap<Integer, IDExpression>();
//         iterator_var = new HashSet<IDExpression>();
        iterator_var = new HashMap<IDExpression, RangeExpression>();
        parameter_map = new HashMap<IDExpression, Integer>();
        iterator_var_map = new HashMap<IDExpression, Integer>();
        total_parameter_ind_map = new HashMap<Integer, IDExpression>();
        total_parameter_spec_map = new HashMap<Integer, Specifier>();
        total_parameter_ind = new HashMap<IDExpression, Integer>();
        global_buffers = new ArrayList<IDExpression>();
        local_buffers = new ArrayList<IDExpression>();
        scalars = new ArrayList<IDExpression>();
    }
    
    void add_parameter(IDExpression expr) {
        parameter.add(expr);
    }
    
    void add_parameter_map(IDExpression expr, Specifier p_spec, int ind, boolean is_prim) {
        if(is_prim)
            parameter_map.put(expr, ind);
        total_parameter_ind_map.put(ind, expr);
        
        total_parameter_spec_map.put(ind, p_spec);
        
        total_parameter_ind.put(expr, ind);
    }
    
    IDExpression get_parameter_idexp(int ind) {
        return total_parameter_ind_map.get(ind);
    }
    
    int get_parameter_id(IDExpression expr) {
        return total_parameter_ind.get(expr);
    }
    
    
    void set_iterator_counter(int value) {
        iterator_counter = 0;
    }
    
//     void add_iterator_variable(IDExpression expr, RangeExpression re) {
    void add_iterator_variable(Map.Entry<IDExpression, RangeExpression> entry) {
        iterator_var.put(entry.getKey(), entry.getValue());
        iterator_var_map.put(entry.getKey(), iterator_counter++);
    }
    
    void add_work_item_func(IDExpression expr, int dim) {
        work_item_func.put(expr, dim);
        work_item_map.put(dim, expr);
    }

    IDExpression get_work_item_func(PolyStart start, int dim) {
        if(work_item_map.containsKey(start.getValue()+dim)) {
            return work_item_map.get(start.getValue()+dim);
        }
        else
            return null;

    }

    void set_non_tdim() {
        String[] non_tdim = {"ss_tx", "ss_ty", "ss_tz"};

        for(int dim=0; dim<3; dim++) {
            if(!work_item_map.containsKey(dim))
                work_item_map.put(dim, new NameID(non_tdim[dim]));
        }
    }

    
    boolean is_iterator_variable(IDExpression expr) {
        return iterator_var.containsKey(expr);
    }

    boolean isGLS(IDExpression id_expr) {
        if(global_buffers.contains(id_expr))
            return true;
        if(local_buffers.contains(id_expr))
            return true;
        if(scalars.contains(id_expr))
            return true;
        return false;
    }
    
    int is_primitive(Expression expr) {
        if(expr instanceof Literal)
            return 64;
        
        if(expr instanceof IDExpression) {
            IDExpression id_expr = (IDExpression) expr;
            if(parameter_map.containsKey(id_expr))      // parameter to the function
                return 65+parameter_map.get(id_expr);
            else if(work_item_func.containsKey(id_expr))        //work item function <64
                return work_item_func.get(id_expr);
            else if(iterator_var.containsKey(id_expr))      // iterator variable
                return 96+iterator_var_map.get(id_expr);
            // else if(id_expr.getName().indexOf("host_")!=-1) {
            // }
                
            String id_name = ((IDExpression)expr).getName();
            if(work_item_variables.indexOf(id_name)!=-1)
                return work_item_variables.indexOf(id_name);
            // for(int i=0; i<work_item_variables.length; i++) {
            //     if(id_name.equals(work_item_variables[i])) {
            //         return i;
            //     }
            // }
        }

        else if(expr instanceof FunctionCall) {
            FunctionCall func_call = (FunctionCall) expr;
            IDExpression func_name = (IDExpression)func_call.getName();
            String str_func_name = func_name.toString();
            
            
            if(work_item_func_names.indexOf(str_func_name)!=-1) {
                int dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
                return 3*work_item_func_names.indexOf(str_func_name)+dim;
            }

            // for(int i=0; i<3; i++) {
            //     if(str_func_name.equals(work_item_func_names[i])) {
            //         int dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
            //         return 3*i+dim;
            //     }
            // }
        }
        return -1;
    }

    Expression expandBasicExpression(Expression tt)    // if expression primitve then it converts the workfunction calls to host variables
    {
        int res = -1;

        if(tt instanceof IDExpression) {
            IDExpression sub_expr = (IDExpression)tt;      
            res = is_primitive(sub_expr);
//                     System.out.println("-- "+ sub_expr.toString()+ " "+res);
        }

        else if(tt instanceof FunctionCall) {
            FunctionCall sub_expr = (FunctionCall)tt;          
            res = is_primitive(sub_expr);
        }

        else if(tt instanceof Literal)
            return tt.clone();

        else return null;

        // System.out.println("Primitive Value: " + res);

        if(res==-1) {
            System.out.println("[ERROR] -1: "+tt.toString());
            return null;
        }
        if(res<64) {
            IDExpression id_iter = TranslationVariableNames.getIteratorVariableHost(res, "", "l");
            return id_iter;
            // sub_expr.swapWith(id_iter);
        }
        else if(res < 96) {
            return tt.clone();
        }

        else if(res >=96) {
            return new NameID("host_"+((IDExpression)tt).getName());
            // return tt.clone();
            // return new NameID(((IDExpression)tt).getName());
            // return new NameID(sub_expr.getName()+"l");
            // sub_expr.swapWith(new NameID(sub_expr.getName()+"l"));
        }   

        return null;
    }

    // set_coefficients
    Expression expandHostExpression(Expression x_expr, Map<IDExpression, Expression> definitions_expand) {
        Expression expr = null;
        if(x_expr instanceof IDExpression) {
            if(definitions_expand.containsKey((IDExpression)x_expr)) {
                expr = definitions_expand.get((IDExpression)x_expr).clone();
            }
            else {
                System.out.println("expandExpression(): Definition not found: " + x_expr);
                expr = x_expr.clone();
            }
        }
        else {
            expr = expandExpression(x_expr, definitions_expand);
        }
        

        if( (expr instanceof FunctionCall) || (expr instanceof IDExpression) ) {
            return expandBasicExpression(expr);
        }
        

        DFIterator<Traversable> expr_iterator = new DFIterator<Traversable>(expr);
        expr_iterator.pruneOn(FunctionCall.class);
        expr_iterator.pruneOn(IDExpression.class);
        
        while(expr_iterator.hasNext()) {
            Traversable tt = expr_iterator.next();
            if ((tt instanceof IDExpression) || (tt instanceof FunctionCall)) {
                Expression iter_expr = (Expression) tt;
                Expression basic_host_expand_expr = expandBasicExpression(iter_expr);

                if(basic_host_expand_expr != null) {
                    iter_expr.swapWith(basic_host_expand_expr);
                }

                else {
                    System.out.println("ISSUE: " + iter_expr);
                    return null;    // not possible with non prime vars
                }
            }  
        }
        
        return expr;
    }

    //expand_expression
    Expression expandExpression(Expression expr, Map<IDExpression, Expression> definitions_expand) {
        List<Traversable> children = expr.getChildren();
        Expression dup_expr = expr.clone();
        
        int child_num = 0;
        for(Traversable child : children) {
            
            if(child instanceof IDExpression) {
                IDExpression child_id = (IDExpression) child;
                if(is_iterator_variable(child_id))  {
//                                 dup_expr.setChild(child_num, new NameID("ss_iter"));
//                             continue;
                }
                if(definitions_expand.containsKey(child_id)) {
                    dup_expr.setChild(child_num, definitions_expand.get(child_id).clone());
                }
            }
            
            else if(!(child instanceof FunctionCall) && !(child instanceof Literal) && (child instanceof Expression)) {
                
                Expression sub_expr = (Expression) child;
//                         System.out.println(sub_expr.toString());
                dup_expr.setChild(child_num, expandExpression(sub_expr, definitions_expand));
            }
            ++child_num;
        }
        return dup_expr;
    }
    
}


class IndexConstInfo {
    // Expression access_expr = null;
    Long const_work_func = 0l;

    Expression const_exprs[] = {null, null, null, null, null, null};


    IndexConstInfo() {
    }

    // IndexInfo(Expression ind_expr) {
        // this.access_expr = ind_expr;
    // }

    boolean add_const_expr_map(int position_work_function, Expression _const) {
        const_work_func = const_work_func | (1l << position_work_function);
        // if(const_exprs[position_work_function] == null ) {

        // }
        const_exprs[position_work_function] = _const;
        return true;
    }

}

class Conditionals {
    List<Expression> all_conditionals;
    List<String> convert_ones = new ArrayList(Arrays.asList("txl", "tyl", "tzl", "bxl", "byl", "bzl", "gxl", "gyl", "gzl"));

    Conditionals() {
        all_conditionals = new ArrayList<Expression>();
    }

    int getIndex(Expression expr) {
        return all_conditionals.indexOf(expr);
    }

    Expression getConditional(int index) {
        if(index < all_conditionals.size())
            return all_conditionals.get(index);
        return null;
    }

    boolean addConditional(Expression expr) {
        boolean new_one = false;

        if(getIndex(expr)==-1) {
            new_one = true;
            all_conditionals.add(expr);
        }

        return new_one;
    }

    boolean addConvertOnes(String c_id_name) {
        convert_ones.add(c_id_name);
        return true;
    }

    Long getConditionalsToStmt(ExpressionStatement expr_stmt, Primitive primitive, Map<IDExpression, Expression> definitions_expand) {
        Traversable cur_trav = expr_stmt;
        // System.out.println("getConditionalsForStmt() : " + expr_stmt);
        Long return_status = 0l;

        while(! (cur_trav instanceof TranslationUnit)) {
            Expression condition = null;

            if(cur_trav instanceof IfStatement) {
                condition = ((IfStatement)cur_trav).getControlExpression();
            }
            else if(cur_trav instanceof ForLoop) {
                condition = ((ForLoop)cur_trav).getCondition();
            }
            else if(cur_trav instanceof WhileLoop) {
                condition = ((WhileLoop)cur_trav).getCondition();
            }

            if(condition != null) {
                condition = primitive.expandHostExpression(condition, definitions_expand);
                if(condition != null) {
                    boolean new_one = addConditional(condition);
                    return_status = return_status | (1l << getIndex(condition));
                    
                    // if(new_one) {
                    //     System.out.println("condition added NEW: " + condition);
                    // }
                    // else
                    //     System.out.println("condition added : " + condition);
                }
            }
            

            cur_trav = cur_trav.getParent();
        }
        return return_status;
    }

    public static int getLastSetBitPos(Long n)
    {
        return (int)(Math.log10(n & ~(n-1)) / Math.log10(2));
    }

    void addIndexing(Expression expr, int ind, Expression modification, boolean do_swap, boolean skipped, int add_to_mod, WorkFuncRangeMap wf_range_map) {
        

        
        Expression m_id = null;
        int count_ids = 0;
        String xdl = "";

        // if(!skipped) {
            // System.out.println("CONDITIONALS PROCESSING: " + expr);
            DFIterator<IDExpression> id_iter = new DFIterator<IDExpression>(expr, IDExpression.class);

            while(id_iter.hasNext()) {
                IDExpression _id = id_iter.next();
                count_ids++;

                if(convert_ones.indexOf(_id.getName()) != -1) {
                    m_id = _id;
                    xdl = _id.getName();

                    // if(do_swap)
                    if(do_swap) {
                        if(!skipped)
                            _id.swapWith(new ArrayAccess(_id.clone(), new IntegerLiteral(ind)));
                    }
                    else
                        System.out.println("Else: "+expr + " : " +_id);
                }
                  
            }
        // }

        // else if(do_swap && (modify_se!=null)) {
        //     // System.out.println("SKIPPED CONDITIONALS PROCESSING: " + expr);
        //     DFIterator<ArrayAccess> id_iter = new DFIterator<ArrayAccess>(modify_se, ArrayAccess.class);

        //     while(id_iter.hasNext()) {
        //         ArrayAccess _id = id_iter.next();
        //         count_ids++;

        //         // System.out.println("ArrayName: " + _id.getArrayName().toString());
        //         if(convert_ones.indexOf(_id.getArrayName().toString()) != -1) {
        //             m_id = _id;

        //             // System.out.println("Else: "+expr + " : " +_id);
        //         }

        //         // System.out.println(count_ids);
                  
        //     }
        // }
        

        // if((count_ids==1) && (modify_se != null) && (m_id!=null)) {
        if((count_ids==1) && (modification!=null) && (m_id!=null)) {
            Expression modify_ref = null;
            Expression ref = null;
            if(m_id instanceof IDExpression) {
                modify_ref = new ArrayAccess(m_id.clone(), new IntegerLiteral(1-ind));
                ref = new ArrayAccess(m_id.clone(), new IntegerLiteral(ind));
            }
            else if(m_id instanceof ArrayAccess)
                modify_ref = m_id;

            if((modify_ref == null) || !(expr.equals(m_id) || expr.equals(ref)) ) {
                System.out.println("NOT EQUAL: " + expr + " " + modify_ref +" " + m_id);
                return;
            }

            IDExpression func_name = null;
            int is_min = 1;
            switch(add_to_mod) {
                case -2:
                    modification = OneArith.minusOne(modification.clone());
                case -1:
                    func_name = new NameID("ss_min");
                    break;
                case 0:
                    is_min = 2;
                    break;
                case 2:
                    modification = OneArith.addOne(modification.clone());
                case 1:
                    func_name = new NameID("ss_max");
                    is_min = 0;
                    break;
            }

            if(func_name != null) {
                wf_range_map.setMap(xdl, 1-ind, modification, is_min);
                modification = new FunctionCall(func_name, modify_ref.clone(), modification.clone());
            }
            else if(is_min == 2) {
                wf_range_map.setMap(xdl, 1-ind, modification, is_min);
            }


            // IRTools.replaceAll(modify_se, modify_ref, modification);
        
        }

    }

    Expression convertConditionBasedOnLOp(Expression expr, boolean do_swap, boolean skipped, WorkFuncRangeMap wf_range_map) { 

        Expression dup_expr = expr.clone();

        if(! (expr instanceof BinaryExpression))
            return dup_expr;

        BinaryExpression bin_expr = (BinaryExpression) dup_expr;
        BinaryOperator bin_op = bin_expr.getOperator();

        Map<BinaryOperator, Integer> addToMod = new HashMap<BinaryOperator, Integer>()
        {{
             put(BinaryOperator.COMPARE_LT, -2);
             put(BinaryOperator.COMPARE_GT, 2);
             put(BinaryOperator.COMPARE_LE, -1);
             put(BinaryOperator.COMPARE_GE, 1);
             put(BinaryOperator.COMPARE_EQ, 0);
        }};

        if( (bin_op == BinaryOperator.COMPARE_LT) || (bin_op == BinaryOperator.COMPARE_LE) ) {

            addIndexing(bin_expr.getLHS(), 0, bin_expr.getRHS(), do_swap, skipped, addToMod.get(bin_op), wf_range_map);
            addIndexing(bin_expr.getRHS(), 1, null, do_swap, skipped, addToMod.get(bin_op), wf_range_map);
        }

        else if( (bin_op == BinaryOperator.COMPARE_GT) || (bin_op == BinaryOperator.COMPARE_GE) ) {

            addIndexing(bin_expr.getLHS(), 1, bin_expr.getRHS(), do_swap, skipped, addToMod.get(bin_op), wf_range_map);
            addIndexing(bin_expr.getRHS(), 0, null, do_swap, skipped, addToMod.get(bin_op), wf_range_map);
        }

        else if ( bin_op == BinaryOperator.COMPARE_EQ) {
            Expression _lhs = bin_expr.getLHS().clone();
            BinaryExpression new_lhs = new BinaryExpression(_lhs, BinaryOperator.COMPARE_LE, bin_expr.getRHS().clone());
            addIndexing(new_lhs.getLHS(), 0, bin_expr.getRHS(), do_swap, skipped, addToMod.get(bin_op), wf_range_map);

            _lhs = bin_expr.getLHS().clone();
            BinaryExpression new_rhs = new BinaryExpression(bin_expr.getRHS().clone(), BinaryOperator.COMPARE_LE, _lhs);
            addIndexing(new_rhs.getRHS(), 1, bin_expr.getRHS(), do_swap, skipped, addToMod.get(bin_op), wf_range_map);

            bin_expr.setOperator(BinaryOperator.LOGICAL_AND);
            bin_expr.setLHS(new_lhs);
            bin_expr.setRHS(new_rhs);
        }

        else if( (bin_op == BinaryOperator.LOGICAL_AND) || (bin_op == BinaryOperator.LOGICAL_OR)) {
            Expression new_lhs = convertConditionBasedOnLOp(bin_expr.getLHS(), do_swap, skipped, wf_range_map);
            Expression new_rhs = convertConditionBasedOnLOp(bin_expr.getRHS(), do_swap, skipped, wf_range_map);

            bin_expr.setLHS(new_lhs);
            bin_expr.setRHS(new_rhs);

        }


        return dup_expr;

    }

    // Pair<IfStatement, CompoundStatement> createIf(long c_cond, int skip)

    Triple<IfStatement, CompoundStatement, Boolean> constructNestedIfStmt(long cur_cond, long prev_cond, CompoundStatement root_one, IfStatement prev_if, WorkFuncRangeMap wf_range_map, KernelRepr _kernel, Primitive k_primitive) {
        // Long prev_start = Long.highestOneBit(prev_index);
        // Long and = prev_cond&cur_cond;
        // Long or = prev_cond|cur_cond;

        if(cur_cond == 0l) {
            return new Triple(null, root_one, false);
        }

        Long xor = prev_cond ^ cur_cond;
        Long skip_bin = (xor& ~(xor-1))-1;
        int skip_actual = Long.bitCount(skip_bin&prev_cond);
        System.out.println("SKIP COUNT: " + skip_actual);

        long skip_cond = skip_bin&prev_cond;
        IfStatement prev_if_copy = prev_if;
        
        IfStatement if_stmt = prev_if;
        CompoundStatement inner_comp = root_one;
        CompoundStatement inner_comp_copy = root_one;
        int skip = skip_actual;

        boolean _first = true;
        boolean atleast_one_new = false;
        
        while(cur_cond!=0l) {
            Long now_over = cur_cond & ~(cur_cond-1);

            if(skip==0) {
                int pos = getLastSetBitPos(now_over);
                CompoundStatement prev_inner_comp = inner_comp;
                inner_comp = new CompoundStatement();

                Expression old_cond = convertConditionBasedOnLOp(getConditional(pos), true, false, wf_range_map);
                Expression new_cond = wf_range_map.modifyToKexpr(_kernel, k_primitive, old_cond);
                // System.out.println("NEW CONDITION = "+ old_cond + "->" + new_cond);
                IfStatement new_if = new IfStatement(new_cond, inner_comp);
                atleast_one_new =true;

                if(prev_inner_comp != null) {
                    prev_inner_comp.addStatement(new_if);
                }
                // if(if_stmt == null)
                if(_first)
                    if_stmt = new_if;
            }
            else {
                inner_comp = (CompoundStatement)prev_if.getThenStatement();
                skip--;

                if(skip!=0) {
                    DFIterator<IfStatement> if_iterator = new DFIterator<IfStatement>(inner_comp, IfStatement.class);
                    while(if_iterator.hasNext()) {
                        prev_if = if_iterator.next();
                    }
                }
            }

            _first = false;
            cur_cond ^= now_over;
        }

        skip = skip_actual;

         while(skip_cond!=0l) {
            // Long now_over = cur_cond_copy & ~(cur_cond_copy-1);
            Long now_over = Long.highestOneBit(skip_cond);

            if(skip>0) {
                int pos = getLastSetBitPos(now_over);
                

                convertConditionBasedOnLOp(getConditional(pos), true, true, wf_range_map);
                
                inner_comp_copy = (CompoundStatement)prev_if_copy.getThenStatement();
                skip--;

                if(skip!=0) {
                    DFIterator<IfStatement> if_iterator = new DFIterator<IfStatement>(inner_comp, IfStatement.class);
                    while(if_iterator.hasNext()) {
                        prev_if_copy = if_iterator.next();
                    }
                }
            }

            skip_cond ^= now_over;
        }

        // System.out.println("_IF: "+ if_stmt);

        return new Triple(if_stmt, inner_comp, atleast_one_new);

    }


}

class IndexToAccessInfo {
    List<Expression> access_map;
    List<Expression> original_access_map;
    List<Long> conditional_map;
    List<Strides> stride_info;
    // List<IndexInfo> access_map;
    Long wmap = 0l;
    Long rmap = 0l;
    // List<IndexConstInfo> const_info;
    // Long const_map = 0;
    // HashMap<Integer, Pair<Long, Expression> > const_exprs[9];

    IndexToAccessInfo() {
        original_access_map = new ArrayList<Expression>();
        access_map = new ArrayList<Expression>();
        // const_info = new ArrayList<IndexConstInfo>();
        conditional_map = new ArrayList<Long>();
        stride_info = new ArrayList<Strides>();

        // for(int i=0; i<9; i++) {
        //     // const_exprs[i] = new ArrayList<Expression>();
        //     const_exprs[i] = new HashMap<Integer, Pair<Long, Expression> >();
        // }
    }

    Pair<Integer, Boolean> add_access_expr(Expression original_expr, Expression expr, int is_write, Long cond_status) {
        int index = access_map.indexOf(expr);
        boolean _new = false;
        if(index==-1) {
            // index = access_map.size();   //CHECK LATER
            _new = true;
            access_map.add(expr);
            original_access_map.add(original_expr);
            // const_info.add(new IndexConstInfo());
            conditional_map.add(cond_status);
            stride_info.add(null);

            index = access_map.indexOf(expr);

        }
        else {
            conditional_map.set(index, conditional_map.get(index) & cond_status);
        }
        
        if((is_write&2) !=0) {
            wmap = wmap | (1l << index);
        }
        if((is_write&1) !=0) {
            rmap = rmap | (1l << index);
        }

        // conditional_map.set(index, conditional_map.get(index)|cond_status);

        System.out.println("Conditionals status = " + Long.toBinaryString(conditional_map.get(index)) );

        return new Pair(index, _new);
        // return _new;
    }

    boolean add_const_expr_map(Expression expr, int position_work_function, Expression const_expr) {
        int index = access_map.indexOf(expr);
        return add_const_expr_map(index, position_work_function, const_expr);
    }

    boolean add_const_expr_map(int index, int position_work_function, Expression const_expr) {
        if(index == -1){
            return false;
        }

        // const_info.get(index).add_const_expr_map(position_work_function, const_expr);
        return true;
    }

    boolean set_stride_info(int index, Strides _stride) {
        if(index >= stride_info.size())
            return false;
        stride_info.set(index, _stride);
        return true;
    }

    Long get_const_expr_map(int index) {
        if(index < conditional_map.size())
            return conditional_map.get(index);
        return null;
    }

    boolean const_work_func_present(int index) {
        // IndexConstInfo ind_const_info = get_const_expr_map(index);
        // if(ind_const_info.const_work_func==0)
        if(conditional_map.get(index)==0l)
            return false;
        return true;
    }

}

class TranslationUnitKernelInfo {
    HashMap<Integer, IndexToAccessInfo> ind_to_access_info;
    Conditionals conditionals_kernel;
    HashMap<Integer, List<Expression>> ind_to_access_map;
    HashMap<Integer, Long> ind_to_access_type_wmap;
    HashMap<Integer, Long> ind_to_access_type_rmap;
    HashMap<Integer, Long> ind_to_access_const_map;
    

    HashMap<IDExpression,Integer> param_ind_map;
    Primitive primitive;
    TranslationUnit tran_unit = null;
    Declaration host_tran_decl_ref = null;
    
    
    boolean params_initialized = false;
    
    TranslationUnitKernelInfo() {
    }
    
    TranslationUnitKernelInfo(HashMap<Integer, IndexToAccessInfo> ind_to_access_info, Conditionals conditionals_kernel , HashMap<IDExpression,Integer> param_ind_map, Primitive primitive, TranslationUnit tran_unit, Declaration host_tran_decl_ref) {
        // this.ind_to_access_map = ind_to_access_map;
        // this.ind_to_access_type_wmap = ind_to_access_type_wmap;
        // this.ind_to_access_type_rmap = ind_to_access_type_rmap;
        this.ind_to_access_info = ind_to_access_info;
        this.conditionals_kernel = conditionals_kernel;
        this.param_ind_map = param_ind_map;
        this.primitive = primitive;
        this.tran_unit = tran_unit;
        this.host_tran_decl_ref = host_tran_decl_ref;
    }
    
    void set_info(HashMap<Integer, IndexToAccessInfo> ind_to_access_info, Conditionals conditionals_kernel, HashMap<IDExpression,Integer> param_ind_map, Primitive primitive) {
        // this.ind_to_access_map = ind_to_access_map;
        // this.ind_to_access_type_wmap = ind_to_access_type_wmap;
        // this.ind_to_access_type_rmap = ind_to_access_type_rmap;
        this.ind_to_access_info = ind_to_access_info;
        this.conditionals_kernel = conditionals_kernel;
        this.param_ind_map = param_ind_map;
        this.primitive = primitive;
    }
    
    void set_translation_unit(TranslationUnit tran_unit) {
        this.tran_unit = tran_unit;
    }
    
    void done_intialization(boolean bool) {
        params_initialized = bool;
    }
    
    boolean is_initialized() {
        return params_initialized;
    }
    
    Declaration getRefDeclaration() {
        return host_tran_decl_ref;
    }
    
    
}

public class OneTBtoDPU extends TransformPass
{

    public OneTBtoDPU(Program program) {
            super(program);
            out = new PrintWriter(System.out);
    }

    public String getPassName() {
            return new String("[OneTBtoDPU]");
    }
    
    List<TranslationUnit> tu_kernel_proc = new ArrayList<TranslationUnit>();
    
    List<Object> add_swap_stmts;
    List<Statement> ref_swap_stmts;
    List<Statement> add_comp_swap_stmts;
    List<Long> swap_dependence;
    List<Boolean> is_swap_next;
    List<Integer> skip_swap_stmts;
    
    List<Object> add_to_stmts;
    List<Statement> ref_to_stmts;
    List<Statement> add_comp_to_stmts;
    List<Long> add_dependence;
    List<Boolean> is_add_next;
    List<Integer> skip_add_stmts;
    
    List<Statement> rm_to_stmts;
    List<Statement> rm_comp_to_stmts;
    
    Primitive primitive;
    
    Map<Symbol, Expression> Definitions;
    Map<IDExpression, Expression> definitions_id;
    Map<IDExpression, Expression> definitions_expand;
            
    Map<Symbol, Long> dep_vector;
    
    HashSet<Expression> primitive_ones;
    HashSet<Symbol> iterator_variables;
    
    //HashSet<Expression> originalExp;
    //HashSet<Symbol> originalSym;
    HashMap<Expression, Symbol> originalES;
    
    Statement cur_stmt;
    CompoundStatement cur_comp_stmt;
    Symbol vector_ids[] = new Symbol[64];
    IDExpression vector_idexps[] = new IDExpression[64];
    //group_ids;
    Symbol group_ids[] = new Symbol[64];
    
//         long partition_idx = (long)(1|(1<<1)|(1<<2));
//         long partition_idx = (long)(1<<1);
//     long partition_idx = (long)((1<<3)|(1<<4)|(1<<5));
    long partition_idx = (long) ((1<<9)-1);
    PrintWriter out;
    Procedure proc = null;
    String proc_name ="";
    
    HashMap<IDExpression, Expression> index_map;

    HashMap<Integer, IndexToAccessInfo> ind_to_access_info;
    HashMap<Integer, List<Expression>> ind_to_access_map;
    // HashMap<Integer, List<boolean>> ind_to_access_type_map;
    HashMap<Integer, Long> ind_to_access_type_wmap;
    HashMap<Integer, Long> ind_to_access_type_rmap;
    Conditionals conditionals;

    HashMap<IDExpression,Integer> param_ind_map;
    
    //static Map<Procedure, TranslationUnitKernelInfo> translation_unit_info;
    
    static Map<String, TranslationUnitKernelInfo> translation_unit_info;
    TranslationUnit tran_unit = null;
    
    Map<Statement, RangeDomain> range_map;
    
    Declaration host_tran_decl_ref = null;
    int BARRIER_COUNT = 0;
    
    
    static {
        //translation_unit_info = new HashMap<Procedure, TranslationUnitKernelInfo>();
        
        translation_unit_info = new HashMap<String, TranslationUnitKernelInfo>();
    }
				
    public void start()
    {
    
        TranslationUnit special_tran_unit = new TranslationUnit("special");
        ss_translation_tu(special_tran_unit);
        
        DFIterator<TranslationUnit> tran_unit_iter = new DFIterator<TranslationUnit>(program, TranslationUnit.class);

// 		TranslationUnit tran_unit = null;
		
		
		
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
            if(!file_extension.equals("cl")) {
                System.out.println("SKIP TranslationUnit");
                continue;
            }
            int fName_start = fileName.lastIndexOf('/');    
            
            
            DepthFirstIterator iter = new DepthFirstIterator(tu);
            proc = null;
            Object obj = null;
            
            while(iter.hasNext()) {
			obj = iter.next();
            
			//only manipulting procedures
            
			if(obj instanceof Procedure) {
			
                proc_name = ((Procedure)obj).getName().getName();
                if(proc_name.equals("get_global_id") || proc_name.equals("get_local_id") || proc_name.equals("get_group_id")|| proc_name.equals("barrier")) {
                    tran_unit.removeChild((Procedure)obj);
                    continue;
                }
                
//                 TranslationUnit new_tu = new TranslationUnit("D_"+proc_name);
//                 new_tu.addDeclaration(((Procedure)obj).clone());
//                 program.addTranslationUnit(new_tu);
                
                BARRIER_COUNT = 0;
            
				add_to_stmts = new ArrayList();
				ref_to_stmts = new ArrayList<Statement>();
				add_comp_to_stmts = new ArrayList<Statement>();
				add_dependence = new ArrayList<Long>();
				is_swap_next = new ArrayList<Boolean>();
				skip_swap_stmts = new ArrayList<Integer>();
				
				add_swap_stmts = new ArrayList();
				ref_swap_stmts = new ArrayList<Statement>();
				add_comp_swap_stmts = new ArrayList<Statement>();
				swap_dependence = new ArrayList<Long>();
				is_add_next = new ArrayList<Boolean>();
				skip_add_stmts = new ArrayList<Integer>();
				
				rm_to_stmts = new ArrayList<Statement>();
				rm_comp_to_stmts = new ArrayList<Statement>();
				
				iterator_variables = new HashSet<Symbol>();
				
				
				//originalExp = new HashSet<Expression>();
				//originalSym = new HashSet<Symbol>();
				originalES = new HashMap<Expression, Symbol>();
				
				Definitions = new HashMap<Symbol, Expression>();
				definitions_id = new HashMap<IDExpression, Expression>();
				definitions_expand = new HashMap<IDExpression, Expression>();
				
				
                dep_vector = new HashMap<Symbol, Long>();
				
				//List<Expression> changedExp = new ArrayList<Expression>();
				
				HashSet<Integer> changed_lines = new HashSet<>();
				proc = (Procedure) obj;
				
				// noting the pointer or array parameters
				Map<Symbol, Expression> pointer_params = new HashMap<Symbol, Expression>();
				// noting non pointer paramaters
				HashSet<Symbol> const_params = new HashSet<Symbol>();
				
				
				index_map = new HashMap<IDExpression, Expression>();
				ind_to_access_map = new HashMap<Integer, List<Expression>>();
                ind_to_access_info = new HashMap<Integer, IndexToAccessInfo>();
                ind_to_access_type_wmap = new HashMap<Integer, Long>();
                ind_to_access_type_rmap = new HashMap<Integer, Long>();
                conditionals = new Conditionals();

				param_ind_map = new HashMap<IDExpression, Integer>();
				
				primitive = new Primitive();

				
                for(int i=0; i<6; i++) {
                    vector_ids[i] = null;
                    group_ids[i] = null;
                    vector_idexps[i] = null;
                }
                    
				//group_ids;
// 				Symbol group_ids[] = new Symbol[3];


				Declarator p_decl = proc.getDeclarator();
				System.out.println("\n--------");
				System.out.println("Current procedure: " + proc.getName().getName());
				p_decl.print(out);
				
				
				
// 				new_tran_unit = new TranslationUnit(proc.getName().getName()+"_new.c");
// 				copyTranslatioUnit(special_tran_unit, new_tran_unit);
// 				Procedure proc_copy = proc.clone();
// 				proc_copy.setName("TEMP");
// 				new_tran_unit.addDeclaration(proc_copy);
// 				proc = proc_copy;
				
                
				// out.write("\n");

				range_map = RangeAnalysis.getRanges(proc);
				
				java.util.List<Declaration> params = proc.getParameters();
				int p_index = 0;
				for(Declaration decl: params) {
					if(decl instanceof VariableDeclaration) {
						//System.out.println(decl.toString());
						VariableDeclaration v_decl = (VariableDeclaration)decl;
                        System.out.println("---DECL--- > "+decl+" "+decl.getClass().getName());
                        Specifier p_spec = Specifier.VOID;
                        
                        List<Specifier> specs = v_decl.getSpecifiers();
                        if(specs.size()>0)
                            p_spec = specs.get(0);
                            

						List<Symbol> sym_list = v_decl.getDeclaredSymbols();
						List<IDExpression> decl_ids = v_decl.getDeclaredIDs();
                        int total = v_decl.getNumDeclarators();
                        Iterator<IDExpression> iter_id = decl_ids.iterator();
                        
						for (Symbol sym: sym_list) {
                            IDExpression id_exp = iter_id.next();
							//System.out.print(sym.getSymbolName() + " : ");
							if (!SymbolTools.isPointer(sym) && !SymbolTools.isArray(sym)) {
								//System.out.println("primitive");
								const_params.add(sym);
								primitive.add_parameter(id_exp);
								primitive.add_parameter_map(id_exp, p_spec, p_index, true);
							}
							else {
								//System.out.println("pointer or array");
								pointer_params.put(sym, new InfExpression(-1));
								param_ind_map.put(id_exp, p_index);
								primitive.add_parameter_map(id_exp, p_spec, p_index, false);
							}	
							++p_index;
						}
						
					}
				}
				out.flush();
				System.out.println("--------");
				
				
				//Iterate procedure body

				CompoundStatement comp_stmt = proc.getBody();
				
				/*out.write("@Range Map@\n");
				for(Map.Entry<Statement, RangeDomain> entry: range_map.entrySet()) {
                    out.write(entry.getKey().toString() + "-> " +entry.getValue().toString() + "\n");
				}*/
				
				String input_0_string = processProcedureDeclarator((ProcedureDeclarator)proc.getDeclarator());

				source_source_translation(comp_stmt, true);

                input_0_string += "if(_bd[PARTITION_DIM_GRID]+t_bd[PARTITION_DIM_GRID]>=GRID_LIMIT) {\n"+
                    "// printf(\"###SKIP### %d %lld %lld\\n\", me(), _bd[PARTITION_DIM_GRID]+t_bd[PARTITION_DIM_GRID], GRID_LIMIT);\n"+
                    "\tint BARRIER_COUNT = "+BARRIER_COUNT+";\n"+
                    "\tfor(int ss_b=0; ss_b<BARRIER_COUNT; ss_b++)\n"+
                    "\t\tbarrier_wait(barrier[t_bd[PARTITION_DIM_GRID]]);\n"+
                    "\treturn;"+
                    "}\n";

                Declaration input_0_decl = createDecls.createCodeAnnotDecl(input_0_string);
				
// 				List<IDExpression> tid_idexp = convertDeclarationStmt(proc, tran_unit);
				convertAccess(comp_stmt, 1);
				
				
// 				for(IDExpression id_exp: tid_idexp) {
//                     IRTools.replaceAll(proc, id_exp, new ArrayAccess(id_exp.clone(), new NameID("tid")));
//                 }

//                 proc.getBody().addDeclarationBefore(proc.getBody().getDeclarations().iterator().next(), createDecls.createCodeAnnotDecl("int tid = me();\n"));
//                 proc.getBody().addDeclaration(input_0_decl);
                
                addDeclarationToBody(proc.getBody(), input_0_decl);
                
				
				
				
				proc.setBody(new CompoundStatement());
// 				proc.setName("none");
                tran_unit.removeChild(proc);
				tran_unit.addDeclaration(createDecls.createProcedure(Specifier.VOID, proc_name, new ArrayList<String>(Arrays.asList("INPUT", "TASKLET_ID_WG")), new ArrayList<String>(Arrays.asList("t_bd")), new ArrayList<String>(Arrays.asList("t_td_se")), comp_stmt));
				
				
				/*out.write("*** BEGIN - Def***\n");
				for(Symbol key: Definitions.keySet()) {
					out.write(key.getSymbolName()+"= ");
					Expression e = Definitions.get(key);
					if(e==null)
						out.write("null");
					else					
						e.print(out);
					out.write("\n");
				}
				out.write("-- END --\n");
				out.flush();*/

                /*out.write("*** BEGIN - Def***\n");
				for(IDExpression key: definitions_expand.keySet()) {
					out.write(key.getName()+"= ");
					Expression e = definitions_expand.get(key);
					if(e==null)
						out.write("null");
					else					
						e.print(out);
					out.write("\n");
				}
				out.write("-- END --\n");
				out.flush();
				
				out.write("@LOCAL_ID@\n");
				for(int i=0; i<6; i++) {
					if(vector_ids[i]==null) {
						out.print("null ");
						continue;
					}
					out.write(i+": "+vector_ids[i].getSymbolName()+" ");
					out.write("\n");
				}
				
				out.write("\n@DEPENDENCES@\n");
				for(Map.Entry<Symbol, Long> entry : dep_vector.entrySet()) {
                    out.write(entry.getKey().getSymbolName()+" "+entry.getValue()+" "+Long.toBinaryString(entry.getValue())+"\n");
				
				}
				out.write("@@@\n");*/
				out.write("@@@INDEXES\n");
				
				for(Map.Entry<Integer, IndexToAccessInfo> entry: ind_to_access_info.entrySet()) {
                    out.write(entry.getKey().toString() + ": ");
                    List<Expression> it = entry.getValue().access_map;
                    if(it==null)
                        continue;
                    int ind = 0;
                    for(Expression expr: it) {
                        out.write(expr + "-" +  ind++ + ", ");
                    }
                    out.write("\n");
				}
				
				out.write("\n@PARAM: INDEX@\n");
				for(Map.Entry<IDExpression, Integer> entry: param_ind_map.entrySet()) {
                    out.write(entry.getKey().toString() +" " + entry.getValue() + "\n");
				}
				
				out.write("\n@ITERATOR VARS@\n");
				
				out.write("\n@PARAMS@\n");
				for(Map.Entry<IDExpression, Integer> entry: primitive.parameter_map.entrySet()) {
                    out.write(entry.getKey() + ": " + entry.getValue()+"\n");
				}
				out.write("@PARAMS@\n\n");
				
// 				Iterator<Symbol> iterator_var_it = iterator_variables.iterator();
// 				while(iterator_var_it.hasNext()) {
//                     Symbol sym = iterator_var_it.next();
//                     out.write(sym.getSymbolName()+"\n");
// 				}

                for(Map.Entry<IDExpression, RangeExpression> entry : primitive.iterator_var.entrySet()) {
                    out.write(entry.getKey() + " " + entry.getValue() + " " + entry.getValue().getLB() +" " + entry.getValue().getUB() + "\n");
                }
				
				
				
				out.flush();
				
				
				ss_translation_tu(tran_unit);
				translation_unit_info.put(proc.getName().getName(), new TranslationUnitKernelInfo(ind_to_access_info, conditionals, param_ind_map, primitive, tran_unit, host_tran_decl_ref));
				tran_unit.setOutputFilename(proc.getName().getName()+".c");

                // translation_unit_info.get(proc.getName().getName()).set_translation_unit(tran_unit);
				
// 				Declaration proc_copy = createDecls.createProcedure(Specifier.VOID, "TEMP_"+proc.getName().getName(), null, null, copyTraversable(proc.getBody()));
// 				new_tran_unit.addDeclaration(proc_copy);
// 				tu_kernel_proc.add(new_tran_unit);
				
                
				//List<Expression, Pair<Integer, String>> change_expression = new ArrayList<Expression, Pair<Integer, String>>();
				
                
                //for(Map.Entry<Expression, Pair<Integer, String>> entry: change_expression.entrySet()) {
                
				
				
                  
                    
				
				
				
				/*Iterator<Integer> lines_iter = changed_lines.iterator();
				while(lines_iter.hasNext()) {
					System.out.println("LINES: "+ lines_iter.next());
				}*/
				
				
				/*Iterator<Expression> orig_iter = originalExp.iterator();
				for(Expression chng:changedExp) {
                    IRTools.replaceAll(cur_comp_stmt, orig_iter.next(), chng);
				}*/
				
				out.flush();
				
			}

        }
        }
        }
//                 for(TranslationUnit tu : tu_kernel_proc) {
//                     program.addTranslationUnit(tu);
//                 }
                System.out.println("HI SAMYUKTHA!!! findArrayIndex is success... :)");
    }
    
    
    void ss_translation_tu(TranslationUnit tu) {
        String code_body;
        Declaration proc_decl;
        
        code_body = "if(INPUT==-1) {\n"+
            "cycles = (uint32_t)1;\n"+
            "return 0;\n"+
            "}\n\n"+
            "if(INPUT==0) {\n"+
            "if(me()==0) {\n"+
            "\tdivide_grid();\n"+
            "\tperfcounter_config(COUNT_CYCLES, true);\n"+
            "++INPUT;\n"+
            "\t//printf(\"%d - NTASKLET_WG = %d %d\\n\", me(), NTASKLET_WG[0], NTASKLET_WG[1]);\n"+        //PRINTF
            "\t//printf(\"%d - NTASK_TASKLET = %d-%d, %d-%d\\n\", me(), NTASK_TASKLET[0][0], NTASK_TASKLET[0][1], NTASK_TASKLET[1][0], NTASK_TASKLET[1][1]);\n"+       //PRINTF
            "}\n"+
            "barrier_wait(&my_barrier);\n\n"+
            "return 0;\n"+
            "}\n"+
            "else if(INPUT==1) {\n"+
            "\tif(me()==0)\n"+
            "\t\tinitialize_BT();\n"+
            "\tbarrier_wait(&my_barrier);\n\n"+
            "}\n"+
            "int TASKLET_ID_WG = me();\n\n"+

            "int t_bd[3] = {_bd[0], _bd[1], _bd[2]};\n"+
            "Point t_td_se[3] = {{0, _bs[0]}, {0, _bs[1]}, {0, _bs[2]}};\n\n"+

            "int category = get_WORK_GROUP_ID(me(), t_bd, &TASKLET_ID_WG);\n\n"+
            "//printf(\"%d - TASKLET_ID_WG = %d\\n\", me(), TASKLET_ID_WG);\n"+
            "//printf(\"%d - t_bd = %d %d %d\\n\", me(), t_bd[0], t_bd[1], t_bd[2]);\n"+

            "get_start_end(TASKLET_ID_WG, category, t_td_se);\n"+
            "//printf(\"%d - t_td_se = %d-%d, %d-%d, %d-%d\\n\", me(), t_td_se[0].start, t_td_se[0].end, t_td_se[1].start, t_td_se[1].end, t_td_se[2].start, t_td_se[2].end);\n "+
            "//T_TD_START[me()] = t_td_se[PARTITION_DIM_WG].start;\n"+
            "//T_TD_END[me()] = t_td_se[PARTITION_DIM_WG].end;\n"+
            "//N_WG_ID[me()] = TASKLET_ID_WG;\n";
            
        code_body += proc_name + "(0, TASKLET_ID_WG, t_bd, t_td_se);\n\n"+
                    "if(me()==0) {\n"+
                    "\tcycles = (uint32_t)perfcounter_get();\n"+
                    "\tprintf(\"cycles = %lld\\n\", cycles);\n"+             //PRINTF
                    "}\nbarrier_wait(&my_barrier);\n\n"+
                    "WORK_DONE[me()] = 1;\n\n"+
                    "return 0;";
                    
        proc_decl = createDecls.createProcedure(Specifier.INT, "main", null, null, null, createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclaration(proc_decl);

        code_body = "printf(\"%s:: \", arr_name);"
                        +"\nfor(int i=0; i<size; i++) {"
                        +    "\nprintf(\"%f \", arr[i]);"
                        +"\n}"
                        +"\nprintf(\"\\n\");";
        List<VariableDeclaration> params = new ArrayList<VariableDeclaration>();
        params.add(new VariableDeclaration(new VariableDeclarator(new ArrayList<Specifier>(Arrays.asList(Specifier.CHAR, PointerSpecifier.UNQUALIFIED)), new NameID("arr_name")) ));
        params.add(new VariableDeclaration(new VariableDeclarator(new ArrayList<Specifier>(Arrays.asList(new UserSpecifier(new NameID("__mram_ptr")), Specifier.FLOAT, PointerSpecifier.UNQUALIFIED)), new NameID("arr")) ));
        params.add(new VariableDeclaration(new VariableDeclarator(new ArrayList<Specifier>(Arrays.asList(Specifier.INT)), new NameID("size")) ));
        ProcedureDeclarator p_declr = new ProcedureDeclarator(new NameID("print_dpu_mat"), params);

        CompoundStatement p_c_st = new CompoundStatement();
        p_c_st.addStatement(createDecls.createCodeAnnotStmt(code_body));
        proc_decl = new Procedure(Specifier.VOID, p_declr, p_c_st);
        tu.addDeclarationFirst(proc_decl);

        code_body = "return start_gd[dim]+get_relative_global_id(dim, WG_ID, _td_dim);\n";

        proc_decl = createDecls.createProcedure(Specifier.INT, "get_global_id", new ArrayList<String>(Arrays.asList("dim", "WG_ID", "_td_dim")), null, null, createDecls.createCodeAnnotStmt(code_body) );
        tu.addDeclarationFirst(proc_decl);

        code_body = "int rem_threads = 0;\n"+
                    // "for(int dim=t_dimension-1; dim>0; dim--) {\n"+
                    "\tif(dim==PARTITION_DIM_GRID)\n"+
                    "\t\t rem_threads += WG_ID*_bs[dim];\n"+
                    // "}\n"+
                    "rem_threads += _td_dim;\n\n"+
                    "\t return rem_threads;\n";
        proc_decl = createDecls.createProcedure(Specifier.INT, "get_relative_global_id", new ArrayList<String>(Arrays.asList("dim", "WG_ID", "_td_dim")), null, null, createDecls.createCodeAnnotStmt(code_body) );
        tu.addDeclarationFirst(proc_decl);


        // code_body = "int rem_threads = 0;\n"+
        //             "for(int dim=t_dimension-1; dim>0; dim--) {\n"+
        //             "\t rem_threads += _td[dim]*_bs[dim-1];\n"+
        //             "}\n"+
        //             "rem_threads += _td[0];\n\n"+
        //             "\t return (TASKLET_ID_WG*B_THREADS)+rem_threads;\n";
        // proc_decl = createDecls.createProcedure(Specifier.INT, "get_global_id", new ArrayList<String>(Arrays.asList("TASKLET_ID_WG")), new ArrayList<String>(Arrays.asList("_td")), null, createDecls.createCodeAnnotStmt(code_body) );
        // tu.addDeclarationFirst(proc_decl);
		
        code_body = "if(PARTITION_DIM_WG == -1)\n"+
                    "\treturn;\n\n"+
                    "int n_TASKLETS_WG = NTASKLET_WG[0] + (1-category);\n\n"+

                    "t_td_se[PARTITION_DIM_WG].start = NTASK_TASKLET[category][0]*TASKLET_ID_WG;\n\n"+
                    "if(TASKLET_ID_WG < NTASK_TASKLET[category][1]) {\n"+
                        "\tt_td_se[PARTITION_DIM_WG].start += TASKLET_ID_WG;\n"+
                        "\tt_td_se[PARTITION_DIM_WG].end = 1;\n"+
                    "}\n"+
                    "else {\n"+
                        "\tt_td_se[PARTITION_DIM_WG].start += NTASK_TASKLET[category][1];\n"+
                        "\tt_td_se[PARTITION_DIM_WG].end = 0;\n"+
                    "}\n"+
                    "t_td_se[PARTITION_DIM_WG].end += t_td_se[PARTITION_DIM_WG].start + NTASK_TASKLET[category][0];";
        proc_decl = createDecls.createProcedure(Specifier.VOID, "get_start_end", new ArrayList<String>(Arrays.asList("TASKLET_ID_WG", "category")),null, new ArrayList<String>(Arrays.asList("t_td_se")), createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclarationFirst(proc_decl);
        
        
        code_body = "if(PARTITION_DIM_GRID==-1)\n"+
                    "\treturn 1;\n\n"+
                    "int TBREAK_GRID = NR_TASKLETS - NTASKLET_WG[1];\n\n"+
                    "if(TASKLET_ID < TBREAK_GRID) {\n"+
                    "\tt_bd[PARTITION_DIM_GRID] = TASKLET_ID/NTASKLET_WG[0];\n"+
                    "\t*TASKLET_ID_WG = TASKLET_ID%NTASKLET_WG[0];\n\n"+
                    "if( t_bd[0] < NTASKLET_WG[1] )\n"+
                    "\treturn 0;\n"+
                    "else\n"+
                    "\treturn 1;\n"+
                    "}\n\n"+
                    "else {\n"+
                    "\tt_bd[PARTITION_DIM_GRID] = TASKLET_ID - TBREAK_GRID;\n"+
                    "\t*TASKLET_ID_WG = NTASKLET_WG[0];\n"+
                    "\treturn 0;\n"+
                    "}";
        proc_decl = createDecls.createProcedure(Specifier.INT, "get_WORK_GROUP_ID", new ArrayList<String>(Arrays.asList("TASKLET_ID")), new ArrayList<String>(Arrays.asList("t_bd", "TASKLET_ID_WG")), null, createDecls.createCodeAnnotStmt(code_body));
        
        tu.addDeclarationFirst(proc_decl);
                    
                
        code_body = "if(PARTITION_DIM_GRID==-1) {\n"+
                    "\tNTASKLET_WG[0] = NR_TASKLETS;\n"+
                    "}\n"+
                    "else {\n"+
                    "\tNTASKLET_WG[0] = NR_TASKLETS/dpu_gs[PARTITION_DIM_GRID];\n"+
                    "\tNTASKLET_WG[1] = NR_TASKLETS%dpu_gs[PARTITION_DIM_GRID];\n"+
                    "}\n\n"+
                    "if(NTASKLET_WG[1]!=0) {\n"+
                    "\tdivide_work_group(0);\n"+
                    "}\n"+
                    "\tdivide_work_group(1);";
        proc_decl = createDecls.createProcedure(Specifier.VOID, "divide_grid", null, null, null, createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclarationFirst(proc_decl);
        
        code_body = "if(PARTITION_DIM_WG==-1)\n"+
                    "\treturn;\n\n"+
                    "int n_TASKS = _bs[PARTITION_DIM_WG];\n"+
                    "int n_TASKLETS_WG = NTASKLET_WG[0] + (1-category);\n\n"+
                    "NTASK_TASKLET[category][0] = n_TASKS/n_TASKLETS_WG;\n"+
                    "NTASK_TASKLET[category][1] = n_TASKS%n_TASKLETS_WG;";
        proc_decl = createDecls.createProcedure(Specifier.VOID, "divide_work_group", new ArrayList<String>(Arrays.asList("category")), null, null, createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclarationFirst(proc_decl);
        
        // code_body = "\tfor(int dim=0; dim<3; dim++) {\n"+
        //             "\tB_THREADS *= _bs[dim];\n"+
        //             "\tif(dim != PARTITION_DIM_WG)\n"+
        //             "\t\tT_THREADS *= _bs[dim];\n"+
        //             "}\n";
        // proc_decl = createDecls.createProcedure(Specifier.VOID, "initialize_BT", null, null, null, createDecls.createCodeAnnotStmt(code_body));
        // tu.addDeclarationFirst(proc_decl);

        code_body = "\tfor(int dim=0; dim<3; dim++) {\n"+
                    "\t start_gd[dim] = _bd[dim] * _bs[dim];\n"+
                    "}\n";
        proc_decl = createDecls.createProcedure(Specifier.VOID, "initialize_BT", null, null, null, createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclarationFirst(proc_decl);
        
        
        code_body = "\treturn (a<b) ?a :b;";
        proc_decl = createDecls.createProcedure(Specifier.INT, "min", new ArrayList<String>(Arrays.asList("a", "b")), null, null, createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclarationFirst(proc_decl);
        
        code_body = "\treturn (a<b) ?b :a;";
        proc_decl = createDecls.createProcedure(Specifier.INT, "max", new ArrayList<String>(Arrays.asList("a", "b")), null, null, createDecls.createCodeAnnotStmt(code_body));
        tu.addDeclarationFirst( proc_decl );
        

        
        
        String code_decl = "#include<stdio.h>\n"+
                    "#include <alloc.h>\n"+
                    "#include <mram.h>\n"+
                    "#include <stddef.h>\n"+
                    "#include <defs.h>\n"+
                    "#include <mutex.h>\n"+
                    "#include <barrier.h>\n"+
                    "#include <perfcounter.h>\n"+
                    "#include <stdint.h>\n\n"+
                    "__host __dma_aligned int64_t N_TASKLETS = NR_TASKLETS;\n"+
                    "#define DPU_MRAM_SIZE 10108864\n"+
                    "MUTEX_INIT(my_mutex);\n"+
                    "BARRIER_INIT(my_barrier, NR_TASKLETS);\n\n"+   //"int N_MULTI_WGS = 1;\n\n"+
                    "__host __dma_aligned int64_t cycles = 0;\n"+
                    "__host __dma_aligned int64_t INPUT = 0;\n"+
                    "__host __dma_aligned int64_t GRID_LIMIT = 1;\n\n"+
                    "//__host __dma_aligned int64_t barrier_count[NR_TASKLETS];\n"+
                    "//__host __dma_aligned int64_t T_TD_START[NR_TASKLETS];\n"+
                    "//__host __dma_aligned int64_t T_TD_END[NR_TASKLETS];\n"+
                    "//__host __dma_aligned int64_t N_WG_ID[NR_TASKLETS];\n"+
                    "__host __dma_aligned int64_t WORK_DONE[NR_TASKLETS];\n\n"+
                    "__host __dma_aligned int64_t MULTI_WGS = 0;\n\n"+
                    "__host __dma_aligned int64_t _bd[] = {0, 0, 0};\n"+
                    "__dma_aligned int64_t start_gd[] = {0, 0, 0};\n"+
                    "__host __dma_aligned int64_t _gs[] = {1, 1, 1};\n"+
                    "__host __dma_aligned int64_t _bs[] = {1, 1, 1};\n\n"+
                    "__host __dma_aligned int64_t dpu_gs[] = {1, 0, 0};\n\n"+
                    "int B_THREADS = 0;\n"+
                    "int T_THREADS = 0;\n\n"+
                    "__host __dma_aligned int64_t PARTITION_DIM_GRID = 1;\n"+
                    "__host __dma_aligned int64_t PARTITION_DIM_WG = 0;\n\n"+
                    "int NTASKLET_WG[2] = {0 , 0};\n"+
                    "int NTASK_TASKLET[2][2] = {{0, 0}, {0, 0}};\n\n"+
                    "typedef struct {\n"+
                    "\tint start;\n"+
                    "\tint end;\n"+
                    "} Point;\n"+
                    "__host __dma_aligned int64_t t_dimension = 2;\n\n";               
        Declaration vars_decl = createDecls.createCodeAnnotDecl(code_decl);
        host_tran_decl_ref = vars_decl; 
        tu.addDeclarationFirst(vars_decl);        
        
    }
    
    void copyTranslatioUnit(TranslationUnit ref, TranslationUnit new_one) {
        DepthFirstIterator<Traversable> df_iterator = new DepthFirstIterator(ref);
        
        while(df_iterator.hasNext()) {
            Traversable tt = df_iterator.next();
            
            if(tt instanceof Declaration) {
                new_one.addDeclaration(((Declaration)tt).clone());
            }
        }
    }
    
    CompoundStatement copyTraversable (Traversable ref) {
        CompoundStatement new_one = new CompoundStatement();
        DepthFirstIterator<Traversable> df_iterator = new DepthFirstIterator(ref);
        
        while(df_iterator.hasNext()) {
            Traversable tt = df_iterator.next();
            
            if(tt instanceof Declaration) {
                new_one.addDeclaration(((Declaration)tt).clone());
            }
            
            else if((tt instanceof Statement) && !(tt instanceof DeclarationStatement)) {
                new_one.addStatement(((Statement)tt).clone());
            }
        }
        return new_one;
    }
    
    
    
    
                
    void addDeclarationToBody(CompoundStatement proc_body, Declaration input_0_decl) {
        Declaration start_decl = createDecls.createCodeAnnotDecl("int tid = me();\nint ss_tx = 0, ss_ty = 0, ss_tz = 0");
        proc_body.addDeclarationBefore(proc_body.getDeclarations().iterator().next(), start_decl);
        proc_body.addDeclaration(input_0_decl);
        
        List<VariableDeclarator> v_decls = new ArrayList<VariableDeclarator>();


        // for(Map.Entry<Integer, IndexToAccessInfo> entry: ind_to_access_info.entrySet()) {
        //     if( (entry.getValue() ==null) || (entry.getValue().access_map.size()<2) )
        //         continue;
                
        //     IDExpression id_exp = primitive.total_parameter_ind_map.get(entry.getKey());
            
        //     for(int i=1; i<entry.getValue().access_map.size(); i++) {
        //         v_decls.add(new VariableDeclarator(new NameID(id_exp.getName()+"_offset"+Integer.toString(i)) ));
        //     }
        // }
        
        // List<Specifier> specs = new ArrayList<Specifier>();
        // specs.add(new UserSpecifier(new NameID("__host")));
        // specs.add(new UserSpecifier(new NameID("__dma_aligned")));
        // // specs.add(Specifier.INT);
        // specs.add(Specifier.LONG);

        // tran_unit.addDeclarationBefore(proc, new VariableDeclaration(specs, v_decls) );
        // proc_body.addDeclarationAfter(start_decl, new VariableDeclaration(Specifier.INT, v_decls) );
        
        
    }
    
    void expressionCorrecter(Procedure proc) {  //It will change the dependent variable accesses to array accesses
        DFIterator<Expression> df_id_exp = new DFIterator(proc, Expression.class);
        int count =0;
        List<Triple<Expression, Integer, IDExpression>> change_expression = new ArrayList<Triple<Expression, Integer, IDExpression>>();
        while(df_id_exp.hasNext()) {
            Traversable tr = df_id_exp.next();
            if(tr instanceof Expression) {
                Expression expr = (Expression)tr;

                List<Traversable> child_expr = expr.getChildren();
                int child_num = 0;
                for(Traversable child: child_expr) {
                    if(child instanceof IDExpression) {
                        IDExpression child_id = (IDExpression)child;
                        for (Map.Entry<Expression,Symbol> entry : originalES.entrySet()) {
                            Expression _exp = entry.getKey();
//                                 for(Expression _exp: originalExp) {
                            IDExpression id_exp = (IDExpression)_exp;
                            if(id_exp.equals(child_id)) {
                                // System.out.println("MATCHED "+ expr.toString()+" "+ child.toString());
//                                         change_expression.add(new Triple(expr, child_num, id_exp.getName()));
                                change_expression.add(new Triple(expr, child_num, id_exp));

                            }
                        }
                    }
                    ++child_num;
                }
                /*Statement st = id_exp.getStatement();
                if(st instanceof DeclarationStatement)
                    continue;
                for(Symbol sym: originalSym) {
                        IRTools.replaceSymbol(id_exp, sym, new ArrayAccess(new NameID(((IDExpression)id_exp).getName()), new NameID("ty")));
                }*/
            }
        }
        
        for(Triple<Expression, Integer, IDExpression> entry : change_expression) {
            Expression expr = entry.getFirst();
            Symbol sym = originalES.get(entry.getThird());
            long local_dep = dep_vector.get(sym);
            ArrayAccess new_arr_acc = createDepArrayAccess(new UnaryExpression( UnaryOperator.DEREFERENCE, entry.getThird().clone()), local_dep&partition_idx);
            // expr.setChild(entry.getSecond(), new_arr_acc);
            // new_arr_acc.setParent(expr);
            IRTools.replaceAll(expr, entry.getThird(), new_arr_acc);
        }
    }

    boolean source_source_translation(Traversable ss_target, boolean initial_check) {
        if(ss_target==null)
            return false;
        DFIterator<Traversable> p_iter = new DFIterator<Traversable>(ss_target, Statement.class);
        Stack<Integer> stack_current = new Stack<Integer>();
        
//         DepthFirstIterator<Traversable> p_iter = new DepthFirstIterator<Traversable>(proc.getBody());

        // pointer_params will tell you mapping each symbol with its rvalue at that time
        

        //about curent environment
        int skip_statements = 0;
        cur_stmt= null;
        Statement changed_stmt = null;
        cur_comp_stmt = null;
        
        //CompoundStatement group_changed_stmts = new CompoundStatement();
        List<Statement> group_changed_stmts = new ArrayList<Statement>();
        long group_dep_vector = 0;
        boolean prev_stmt_independent = true;
        int line_num = 0;
        int minus=0;
        int minus_index = -1;
        while(p_iter.hasNext()) {
            Traversable o = p_iter.next();
//             out.write(line_num+" "+o.getClass().getName()+"\n");
            ++line_num;
            if(o instanceof Statement) {
                
                boolean cur_stmt_independent = true;
                cur_stmt = (Statement) o;
                
                //AnnotationStatement, BreakStatement, Case, CompoundStatement, ContinueStatement, DeclarationStatement, Default, DoLoop, ExpressionStatement, ForLoop, GotoStatement, IfStatement, Label, NullStatement, ReturnStatement, SwitchStatement, WhileLoop
                
                //System.out.println(":: "+ cur_stmt.getClass().getName()+" " + temp.getParent().getClass().getName()+ " " + cur_stmt.where() +"\n");	
                
                if(cur_stmt.getParent() instanceof CompoundStatement){
                    cur_comp_stmt = (CompoundStatement)cur_stmt.getParent();
                }
                if(initial_check && (o instanceof ExpressionStatement)) {
                    Expression expr = (Expression)((ExpressionStatement)o).getExpression();
                    findArrayAccessPattern(expr);
                    if(expr instanceof AssignmentExpression) {
                        AssignmentExpression ass_expr = (AssignmentExpression)expr;
                        Expression lhs = ass_expr.getLHS();;	
                        if(lhs instanceof IDExpression) {
                            Set<Symbol> lhs_syms = SymbolTools.getAccessedSymbols(lhs);
                            for(Symbol lhs_sym: lhs_syms) {					
                                Definitions.put(lhs_sym, ass_expr.getRHS());
                                definitions_id.put((IDExpression)lhs, ass_expr.getRHS());
//                                 if(! primitive.is_iterator_variable((IDExpression)lhs)) {
                                if(!iterator_variables.contains(lhs_sym)) {
                                    definitions_expand.put((IDExpression)lhs, expand_expression(ass_expr.getRHS()));
                                }
                                else {
                                    System.out.println("12-05-2021: "+lhs.toString());
                                }
                                long local_dep = getDependenceVector(ass_expr.getRHS(), dep_vector, proc);
                                if(dep_vector.containsKey(lhs_sym))
                                    dep_vector.put(lhs_sym, dep_vector.get(lhs_sym)|local_dep);
                                else
                                    dep_vector.put(lhs_sym, local_dep);
                                break;
                            }
                        }
                        
                        /*List<Traversable> child_expr = ass_expr.getChildren();
                        for(Traversable tt: child_expr) {
                            // out.write(tt.getClass().getName()+" ");
                        }*/
    
                    }
                }
                
                
                if(o instanceof IfStatement) {
                    
                    IfStatement if_stmt = (IfStatement)o;
                    Expression cond = if_stmt.getControlExpression();
                    long local_dep = getDependenceVector(cond, dep_vector, proc);
// 							|getDependenceVector(if_stmt.getThenStatement(), dep_vector, proc)|getDependenceVector(if_stmt.getElseStatement(), dep_vector, proc);
                    if(local_dep!=0) {
                        swap_with_statement(if_stmt, if_stmt, cur_comp_stmt, local_dep, prev_stmt_independent);
                        out.write("IF: "+Long.toBinaryString(local_dep)+"\n");
                        
                        skip_statements += numSkipStatements(if_stmt);
                        stack_current.push(ref_swap_stmts.size()-1);
                        minus = 0;
                        line_num += skip_statements;
                        
                        cur_stmt_independent = false;
                        continue;
                    }
                    out.flush();
                    
                }
                
                else if(o instanceof ForLoop || o instanceof WhileLoop || o instanceof DoLoop) {
                    
                    Statement st=null;
                    Expression cond = null;
                    if(o instanceof ForLoop) {
                        cond = ((ForLoop) o).getCondition();
                        primitive.add_iterator_variable(findLoopIndex((ForLoop) o));
                    }
                    else if(o instanceof WhileLoop)
                        cond = ((WhileLoop) o).getCondition();
                    else
                        cond = ((DoLoop) o).getCondition();
                        
                    long local_dep = getDependenceVector(cond, dep_vector, proc);
                    if(local_dep!=0) {
                        skip_statements += numSkipStatements((Statement)o);
                        line_num += skip_statements;
                        
                        swap_with_statement(cur_stmt, cur_stmt, cur_comp_stmt, local_dep, prev_stmt_independent);
                        stack_current.push(ref_swap_stmts.size()-1);
                        minus = 0;
//                         minus_index = ref_swap_stmts.size()-1;
                        cur_stmt_independent = false;
                        continue;
                    }
                    
                }

                else if(o instanceof DeclarationStatement) {
                
//                         cur_stmt_independent = !(getDefinitionsAndDependence(proc, o, prev_stmt_independent));
                        cur_stmt_independent = !(getDefinitionsAndDependence(proc, o, true));
                    
                }
                
                else if(o instanceof ExpressionStatement) {
                            
                    long local_dep = getDependenceVector(cur_stmt, dep_vector, proc);

                    //if(!changed_lines.contains(cur_stmt.where()) && (local_dep & partition_idx)!=0) {
                    if(local_dep!=0) {
                        
                        //changed_stmt = cur_stmt.clone();
                        changed_stmt = cur_stmt;
                        
                        swap_with_statement(cur_stmt, cur_stmt, cur_comp_stmt, local_dep, prev_stmt_independent);
                        
                        //add_to_stmts.add(for_loop);
                        
                        //changed_lines.add(cur_stmt.where());
                        //comp_stmt.removeStatement(cur_stmt);		
                    }
                    
                    out.flush();
                }
                if(skip_statements>0) {
                    if(!cur_stmt_independent) {
                        ++minus;
                    }
                    if((skip_statements==1) && (stack_current.size()>0)) {
                        int temp = stack_current.pop();
                        skip_swap_stmts.set(temp, minus);
                        
                    }
                    skip_statements--;
                }
                prev_stmt_independent = cur_stmt_independent;
            }
            
        }

        
        expressionCorrecter(proc);
        
        
        add_statements_IR();

        swap_statements_IR();

        primitive.set_non_tdim();

        // for(Map.Entry<Integer, IndexToAccessInfo> entry: ind_to_access_info.entrySet()) {
        //     for(Expression expr: entry.getValue().original_access_map) {
        //         List<Expression> find_expr = IRTools.findExpressions(proc, expr);
        //         System.out.println(expr+": "+find_expr.size());
        //         for(Expression _expr: find_expr) {
        //             if(_expr == expr)
        //                 System.out.print(" Y ");
        //             else
        //                 System.out.print(" N ");
        //         }
        //     }
        //     System.out.println("----------");
        // }


        addRelativeScalars(proc, primitive);
         
        convertWorkFunctionCalls(proc, primitive, false);


        return true;
    }
    
//     void convertWorkFunctionCalls(Procedure proc) {
// //         DFIterator<FunctionCall> df_func_iterator = new DFIterator<FunctionCall>(proc, FunctionCall.class);
    
// //         DFIterator<Traversable> df_func_iterator = new DFIterator<Traversable>(proc, Traversable.class);
//         DepthFirstIterator df_func_iterator = new DepthFirstIterator(proc);
//         CompoundStatement c_st = null;
        
//         while(df_func_iterator.hasNext()) {
//             Object tt = df_func_iterator.next();
//             if(tt instanceof CompoundStatement) {
//                 c_st = (CompoundStatement)tt;
//             }

//             if(tt instanceof IDExpression) {

//             }
            
//             if(tt instanceof FunctionCall) {
//                 FunctionCall func_call = (FunctionCall)tt;
//                 Expression call_name = func_call.getName();
                
//                 int offset = -1, dim = 0;
//                 Statement st = null;

//                 Statement stmt = func_call.getStatement();
//                 //out.write("CALL_NAME: "+call_name.toString()+" ");
//                 if(call_name.toString().equals("get_local_id")) {
//                     offset=0;
//                 }
//                 else if(call_name.toString().equals("get_global_id")) {
//                     offset = 3;
//                 }
//                 else if(call_name.toString().equals("get_group_id")) {
//                     offset = 6;
//                 }  
//                 else if(call_name.toString().equals("barrier")) {
                    
//                     Statement new_stmt = createDecls.createCodeAnnotStmt("barrier_wait(barrier[t_bd[PARTITION_DIM_GRID]]);\n++barrier_count[tid];\n");
//                     stmt.swapWith(new_stmt);
// //                     c_st.addStatementBefore(new_stmt, createDecls.createCodeAnnotStmt("++barrier_count[tid];"));
//                 }
//                 else if(call_name.toString().equals("printf")) {
//                     stmt.swapWith(new NullStatement());
//                 }
                
//                 if(offset!=-1) {
//                     dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
                    
//                 }
                
//                 switch(offset) {
//                     case 0:
// //                         func_call.swapWith(new ArrayAccess(getIteratorVariable(offset+dim).clone(), func_call.getArgument(0).clone()));
// //                         st = func_call.getStatement();
// //                         if(st instanceof DeclarationStatement) {
//                             func_call.swapWith(new IntegerLiteral(0));
// //                         }
//                         break;
//                     case 3:
//                             if((stmt instanceof IfStatement) || (stmt instanceof WhileLoop)) {
//                                 func_call.swapWith(new FunctionCall(new NameID("get_relative_global_id"), new NameID("TASKLET_ID_WG"), new NameID("_td")));
//                             }
//                             else
//                                 func_call.swapWith(new FunctionCall(new NameID("get_global_id"), new NameID("TASKLET_ID_WG"), new NameID("_td")));
                                
//                                 List<Expression> wid_exprs= new ArrayList<Expression>();
// //                                 long local_dep = 1l | (1l<<1);
//                                 for(int t_dim=0; t_dim<3; t_dim++) {
//                                     wid_exprs.add(getIteratorVariable(t_dim));
//                                 }
// //                                 while(local_dep!=0) {
// //                                     wid_exprs.add(getIteratorVariable(local_dep& ~(local_dep-1)));
// //                                     local_dep = local_dep&(local_dep-1);
// //                                 }
//                                 Declaration td_decl = createArrayVariableDeclaration(new NameID("_td"), wid_exprs);
//                                 c_st.addDeclaration(td_decl);

// //                             func_call.swapWith(new BinaryExpression(new ArrayAccess(new NameID("t_bd"), func_call.getArgument(0).clone()), BinaryOperator.MULTIPLY, getIteratorLimit(dim)));
// //                         }
//                         break;
                        
//                     case 6:
//                         func_call.swapWith(new ArrayAccess(new NameID("t_bd"), func_call.getArgument(0).clone()));
//                 }
//             }
//         }
				
//     }

    boolean isExprBelongsToCondition(Statement stmt, Expression expr) {
        if((stmt instanceof IfStatement) || (stmt instanceof WhileLoop))
            return true;
        if(stmt instanceof ForLoop) {
            ForLoop for_loop = (ForLoop)stmt;
            if(IRTools.containsExpression(for_loop.getCondition(), expr))
                return true;
        }

        return false;
    }
    void convertWorkFunctionCalls(Traversable trav, Primitive primitive, boolean change_irrespective) {

        DepthFirstIterator df_func_iterator = new DepthFirstIterator(trav);
        CompoundStatement c_st = null;
        ArrayAccess cur_array_access = null;
        IDExpression cur_array_name = null;
        while(df_func_iterator.hasNext()) {
            Object tt = df_func_iterator.next();
            if(tt instanceof CompoundStatement) {
                c_st = (CompoundStatement)tt;
            }

            if(tt instanceof ArrayAccess) {
                cur_array_access = (ArrayAccess)tt;

                // DFIterator<IDExpression> df_array_iter = new DFIterator<IDExpression>(cur_array_access.getArrayName(), IDExpression.class);
                // if(df_array_iter.hasNext()) {
                    if(!(cur_array_access.getArrayName() instanceof IDExpression))
                        continue;
                    IDExpression temp_array_name = (IDExpression)cur_array_access.getArrayName();
                    // IDExpression temp_array_name = df_array_iter.next();
                    // if(primitive.isGLS(temp_array_name)) {
                    if(primitive.global_buffers.indexOf(temp_array_name)!=-1) {
                        cur_array_name = temp_array_name;


                        // if(primitive.total_parameter_ind.containsKey(cur_array_name)) {
                        //     // System.out.println("Before: ");
                        //     for(Expression check_expr : ind_to_access_info.get(primitive.total_parameter_ind.get(cur_array_name)).original_access_map) {
                        //         System.out.println(check_expr);
                        //         if(check_expr.equals(cur_array_access.getIndex(0))) {
                        //             if(check_expr == cur_array_access.getIndex(0))
                        //                 System.out.print(" YES ");
                        //             else
                        //                 System.out.print(" NO ");
                        //         }
                        //     }
                        //     System.out.println(": "+cur_array_access.getIndex(0));
                        // }
                        
                        
                        DFIterator<Expression> df_ind_iter = new DFIterator<Expression>(cur_array_access.getIndex(0), Expression.class);
                        IDExpression modify_array_name = null;

                        while(df_ind_iter.hasNext()) {
                            Expression _expr = df_ind_iter.next();

                            if(_expr instanceof ArrayAccess) {
                                modify_array_name = IRTools.getExpressionsOfType(_expr, IDExpression.class).get(0);
                                if(primitive.isGLS(modify_array_name)){
                                    modify_array_name = null;
                                }
                            }

                            else if(_expr instanceof IDExpression) {
                                IDExpression cur_id = (IDExpression)_expr;
                                if(modify_array_name==null || (!cur_id.equals(modify_array_name)))
                                    continue;
                                // if((cur_array_name==null)||change_irrespective)
                                    // continue;

                                int prim_ind = primitive.is_primitive(cur_id);

                                if( !(cur_id.equals(cur_array_name)) && IRTools.isAncestorOf(cur_array_access, cur_id)) {
                                    if((prim_ind >2) && (prim_ind < 64))
                                        cur_id.swapWith(new NameID(cur_id.getName() + "_relative"));

                                    else if(prim_ind==-1) {
                                        if(dep_vector.containsKey(cur_id)) {
                                            long local_dep = dep_vector.get(cur_id);
                                            if(local_dep!=0l)
                                                cur_id.swapWith(new NameID(cur_id.getName() + "_relative"));
                                        }
                                    }
                                }
                            }
                        }
                        

                        // if(primitive.total_parameter_ind.containsKey(cur_array_name)) {
                            // System.out.println("After: ");
                            // for(Expression check_expr : ind_to_access_info.get(primitive.total_parameter_ind.get(cur_array_name)).original_access_map) {
                            //     System.out.println(check_expr);
                            // }
                            // System.out.println(": "+cur_array_access.getIndex(0));
                        // }
                    }
                // }

            }
            // if(tt instanceof IDExpression) {
            //     IDExpression cur_id = (IDExpression)tt;
            //     if((cur_array_name==null)||change_irrespective)
            //         continue;

            //     int prim_ind = primitive.is_primitive(cur_id);

            //     if( !(cur_id.equals(cur_array_name)) && IRTools.isAncestorOf(cur_array_access, cur_id)) {
            //         if((prim_ind >2) && (prim_ind < 64)) {
            //             cur_id.swapWith(new NameID(cur_id.getName() + "_relative"));
            //         }
            //         // else if(prim_ind==-1) {
            //         //     if(dep_vector.containsKey(cur_id)) {
            //         //         long local_dep = dep_vector.get(cur_id);
            //         //         if(local_dep!=0l) {
            //         //             cur_id.swapWith(new NameID(cur_id.getName() + "_relative"));
            //         //         }
            //         //     }
            //         // }
            //     }
            // }
            
            if(tt instanceof FunctionCall) {
                FunctionCall func_call = (FunctionCall)tt;
                Expression call_name = func_call.getName();
                
                int offset = -1, dim = 0;
                Statement st = null;

                Statement stmt = func_call.getStatement();
                //out.write("CALL_NAME: "+call_name.toString()+" ");
                if(call_name.toString().equals("get_local_id")) {
                    offset=0;
                }
                else if(call_name.toString().equals("get_global_id")) {
                    offset = 3;
                }
                else if(call_name.toString().equals("get_group_id")) {
                    offset = 6;
                }  
                else if(call_name.toString().equals("barrier")) {
                    BARRIER_COUNT += 1;
                    Statement new_stmt = createDecls.createCodeAnnotStmt("barrier_wait(barrier[t_bd[PARTITION_DIM_GRID]]);\n//++barrier_count[tid];\n");
                    stmt.swapWith(new_stmt);
//                     c_st.addStatementBefore(new_stmt, createDecls.createCodeAnnotStmt("++barrier_count[tid];"));
                }
                else if(call_name.toString().equals("printf")) {
                    stmt.swapWith(new NullStatement());
                }
                
                if(offset!=-1) {
                    dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
                    
                }
                
                switch(offset) {
                    case 0:
                        // if(IRTools.isAncestorOf(cur_array_access, func_call)) {
                            IDExpression mod_id = primitive.get_work_item_func(PolyStart.T_START, dim);
                            func_call.swapWith(mod_id.clone());
                        // }
                            // func_call.setFunction(new NameID(func_call.getName().toString()+"_relative"));
                            // func_call.swapWith(new IntegerLiteral(0));
                        break;
                    case 3:
                        if(IRTools.isAncestorOf(cur_array_access, func_call))
                            func_call.swapWith(new FunctionCall(new NameID("get_relative_global_id"), new IntegerLiteral(dim), new ArrayAccess(new NameID("t_bd"), new IntegerLiteral(dim)), primitive.get_work_item_func(PolyStart.T_START,dim).clone()));
                        else
                            func_call.swapWith(new FunctionCall(new NameID("get_global_id"), new IntegerLiteral(dim), new ArrayAccess(new NameID("t_bd"), new IntegerLiteral(dim)), primitive.get_work_item_func(PolyStart.T_START, dim).clone()));

                            // func_call.swapWith(new FunctionCall(new NameID("get_global_id"), new IntegerLiteral(dim), new NameID("TASKLET_ID_WG"), new NameID("_td"))); 
                        
                        List<Expression> wid_exprs= new ArrayList<Expression>();
//                                 long local_dep = 1l | (1l<<1);
                        for(int t_dim=0; t_dim<3; t_dim++) {
                            wid_exprs.add(getIteratorVariable(t_dim));
                        }

                        // Declaration td_decl = createArrayVariableDeclaration(new NameID("_td"), wid_exprs);
                        // c_st.addDeclaration(td_dSecl);

                        break;
                        
                    case 6:
                        if(IRTools.isAncestorOf(cur_array_access, func_call))
                            func_call.swapWith(new ArrayAccess(new NameID("t_bd"), func_call.getArgument(0).clone()));
                        else {
                            Expression _lhs = new ArrayAccess(new NameID("_bd"), func_call.getArgument(0).clone());
                            Expression _rhs = new ArrayAccess(new NameID("t_bd"), func_call.getArgument(0).clone());
                            func_call.swapWith(new BinaryExpression(_lhs, BinaryOperator.ADD, _rhs));
                        }
                
                }
            }
        }
                
    }
    
    
    
    List<IDExpression> convertDeclarationStmt(Procedure tt, TranslationUnit tran_unit) {
    
        List<IDExpression> tid_idexp = new ArrayList<IDExpression>();
    
        DFIterator<DeclarationStatement> df_iterator = new DFIterator<DeclarationStatement>(tt, DeclarationStatement.class);
        
        while(df_iterator.hasNext()) {
            DeclarationStatement decl_stmt = df_iterator.next();
            
//             if(! (decl_stmt.getDeclaration() instanceof VariableDeclaration) || !(decl_stmt.getParent().equals(proc.getBody())) )
            if(! (decl_stmt.getDeclaration() instanceof VariableDeclaration) || (decl_stmt.getParent() instanceof ForLoop))
                continue;
                
            VariableDeclaration decl = (VariableDeclaration) decl_stmt.getDeclaration();
            System.out.println("Decl: " + decl + " " + decl_stmt.getParent().getClass().getName() + " " + getDependenceVector(decl, dep_vector, tt));
            
//             tid_idexp.addAll(decl.getDeclaredIDs());
            
            int num_declarators = decl.getNumDeclarators();
            
//             List<ExpressionStatement> declaration_init = new ArrayList<ExpressionStatement>();
            CompoundStatement declaration_init = new CompoundStatement();
            
            for(int i=0; i<num_declarators; i++) {
                VariableDeclarator v_declr = (VariableDeclarator) decl.getDeclarator(i);
                
//                 long local_dep = getDependenceVector(v_declr, dep_vector, tt);
                
//                 if((local_dep>0) && (local_dep<3))
//                     continue;
//                 boolean not_add = false;
//                 for(int ind=0; ind<3; ind++) {
//                     if((vector_idexps[ind]!=null) && (vector_idexps[ind].equals(v_declr.getID())) )
//                     {
//                         not_add = true;
//                         break;
//                     }
//                 }
//                 if(not_add)
//                     continue;
//                 
//                 System.out.println("::" + v_declr.getID());
//                     
//                 tid_idexp.add(v_declr.getID());
                
                List<Specifier> t_specs = v_declr.getArraySpecifiers();
                
                Expression add_trail_spec = new NameID("N_TASKLETS");
                if(t_specs.size() != 0) {
//                     List<Expression> dim_exprs = new ArrayList<Expression>();
                    ArraySpecifier arr_spec = (ArraySpecifier)t_specs.get(0);
                    
                    for(int j=0; j<arr_spec.getNumDimensions(); j++) {
                        Expression temp = arr_spec.getDimension(j);
                        arr_spec.setDimension(j, add_trail_spec);
                        add_trail_spec = temp;
                    }
                    
//                     t_specs.set(0, new ArraySpec ifier(dim_exprs));
//                     ((ArraySpecifier)t_specs.get(0)).setDimension(0, new NameID("N_TASKLETS"));
                }
                else {
                    tid_idexp.add(v_declr.getID());
                }
                
                
                v_declr.addTrailingSpecifier(new ArraySpecifier(add_trail_spec));
                Initializer init = v_declr.getInitializer();
                if(init!=null) {
                        List init_child = init.getChildren();
                        if(init_child.size()==1)
                            declaration_init.addStatement(new ExpressionStatement(new AssignmentExpression(v_declr.getID().clone(), AssignmentOperator.NORMAL, ((Expression)(init_child.get(0))).clone())));
                        v_declr.setInitializer(null);
//                     v_declr.swapWith(new AssignmentExpression());
                }
            }
            
    //         Declaration decl_copy = decl.clone();
//             List<IDExpression> decl_ids = decl.getDeclaredIDs();
//             List<VariableDeclarator> declr_ids = new ArrayList<VariableDeclarator>();
            
//             for(IDExpression decl_id : decl_ids) {
//                 declr_ids.add(new VariableDeclarator(decl_id.clone()));
//             }

            decl_stmt.swapWith(declaration_init);
            tran_unit.addDeclarationBefore(proc, decl);
            
//                 proc.getBody().addStatement();
//             tran_unit.addDeclarationBefore(proc, new VariableDeclaration(decl.getSpecifiers(), declr_ids));
                
        }
        
        
        
        System.out.println("Declaration converted!!!");
        return tid_idexp;
        
    }

    void addRelativeScalars(Procedure proc, Primitive primitive) {
        DFIterator<Statement> dfs_iter = new DFIterator<Statement> (proc, Statement.class);
        // dfs_iter.pruneOn(Statement.class);
        CompoundStatement c_st = null;
        int skip = 0;
        while(dfs_iter.hasNext()) {
            Statement o_stmt = dfs_iter.next();

            if(o_stmt instanceof CompoundStatement)
                c_st = (CompoundStatement)o_stmt;
            else if(o_stmt.getParent() instanceof CompoundStatement)
                c_st = (CompoundStatement)(o_stmt.getParent());

            // System.out.println("addRelativeScalars(): "+o_stmt);
            if(skip>0) {
                // System.out.println("SKIPPED: " +o_stmt + " " + o_stmt.getClass().getName());

                skip--;
                continue;
            }

            if(o_stmt instanceof DeclarationStatement)  {
                Statement stmt = o_stmt;

                Declaration decl = ((DeclarationStatement)stmt).getDeclaration();

                // System.out.println("addRelativeScalars-decl()0: "+decl + " " + decl.getClass().getName());

                if(! (decl instanceof VariableDeclaration) )
                    continue;
                VariableDeclaration v_decl = (VariableDeclaration)decl;
                
                if(v_decl.getNumDeclarators()>1) {
                    System.out.println("[ERROR] more declarators in one declaration");
                    continue;
                }
                Declarator declr = v_decl.getDeclarator(0); 
                // System.out.println("addRelativeScalars-decl()1: "+o_stmt);
                if(declr instanceof NestedDeclarator) {
                    // System.out.println("addRelativeScalars()-decl2: "+o_stmt);
                    VariableDeclarator v_declr = (VariableDeclarator)((NestedDeclarator)declr).getDeclarator();
                    if( !(primitive.isGLS(v_declr.getID()) ) ){
                        // System.out.println("addRelativeScalars-decl()3: "+o_stmt);
                        String code = stmt.toString();

                        // System.out.println("Before Replace: " + code);
                        code = code.replaceFirst(v_declr.getID().getName(), v_declr.getID().getName()+"_relative");
                        // System.out.println("After Replacing  "+ v_declr.getID().getName()+": " + code);
                        // v_declr.setName(v_declr.getID().getName()+"_relative");
                        // proc.getBody().addDeclarationAfter(((DeclarationStatement)o_stmt).getDeclaration(), ((DeclarationStatement)stmt).getDeclaration());
                        c_st.addDeclarationAfter(((DeclarationStatement)o_stmt).getDeclaration(), createDecls.createCodeAnnotDecl(code));
                        // skip ++;
                    }
                }
            }

            else if(o_stmt instanceof ExpressionStatement) {
                Expression expr = ((ExpressionStatement)o_stmt).getExpression();

                if(expr instanceof AssignmentExpression) {
                    AssignmentExpression assign_expr = ((AssignmentExpression)expr).clone();

                    if(assign_expr.getLHS() instanceof ArrayAccess) {
                        ArrayAccess arr_acc = (ArrayAccess)assign_expr.getLHS();
                        IDExpression cur_array_name = null;

                        DFIterator<IDExpression> df_array_iter = new DFIterator<IDExpression>(arr_acc.getArrayName(), IDExpression.class);
                        if(df_array_iter.hasNext())
                            cur_array_name = df_array_iter.next();
                        // System.out.println("addRelativeScalars()2: "+cur_array_name);
                        if( !(primitive.isGLS(cur_array_name)) ){
                            

                            // int res = removeRedundants(assign_expr.getRHS());
                            // if(res==-1)
                            //     System.out.println("[ERROR] check index access pattern "+ assign_expr.getRHS());
                            
                            

                            cur_array_name.swapWith(new NameID(cur_array_name.getName()+"_relative"));
                            // System.out.println("addRelativeScalars()3: "+o_stmt);
                            // System.out.println("addRelativeScalars()4: "+assign_expr);
                            convertWorkFunctionCalls(assign_expr.getRHS(), primitive, true);
                            c_st.addStatementAfter(o_stmt, new ExpressionStatement(assign_expr));
                            // skip++;
                        }

                    }
                }
            }

        }
    }

    int removeRedundants(Expression expr, boolean replace_0) {     //remove the constant terms
        if(expr == null)
            return 0;
        if(expr instanceof Literal)
            return 0;
        else if(expr instanceof IDExpression) {
            if(primitive.scalars.contains((IDExpression)expr))
                return 0;
            return 1;
        }
        else if(expr instanceof FunctionCall) {
            return 1;
        }
        else if(expr instanceof ArrayAccess) {
            return 1;
        }
        else if(expr instanceof UnaryExpression) {
            UnaryOperator u_op = ((UnaryExpression)expr).getOperator();
            if((u_op == UnaryOperator.BITWISE_COMPLEMENT)
                || (u_op == UnaryOperator.LOGICAL_NEGATION)
                || (u_op == UnaryOperator.MINUS)
                || (u_op == UnaryOperator.PLUS))
                return removeRedundants(((UnaryExpression)expr).getExpression(), replace_0);
            return -1;
        }
        if(!(expr instanceof BinaryExpression))
            return -1;

        BinaryExpression bin_expr = (BinaryExpression)expr;

        BinaryOperator bin_op = bin_expr.getOperator();
        Expression bin_lhs = bin_expr.getLHS();
        Expression bin_rhs = bin_expr.getRHS();

        if(replace_0) {
            if(bin_op != BinaryOperator.ADD && bin_op != BinaryOperator.SUBTRACT)
                replace_0 = false;
        }
        
        int _lhs = removeRedundants(bin_lhs, replace_0);
        int _rhs = removeRedundants(bin_rhs, replace_0);

        if(_lhs==-1 || _rhs==-1)
            return -1;

        if( bin_op== BinaryOperator.ADD || bin_op==BinaryOperator.SUBTRACT) {
            if(_lhs==0 && replace_0) {
                bin_lhs.swapWith(new IntegerLiteral(0));
            }
            if(_rhs==0 && replace_0)
                bin_rhs.swapWith(new IntegerLiteral(0));

            return Math.min(_lhs+_rhs, 1);
        }

        if(bin_op == BinaryOperator.MULTIPLY)
            if(_lhs+_rhs > 1)
                return -1;
        return _lhs+_rhs;

        // if(_lhs==1 && _rhs==1)
            // return 1;
        // return 0;
    }

    
    
    void convertAccess(Traversable tt, int local) {
        DFIterator<Expression> df_ex_iterator = new DFIterator<Expression>(tt, Expression.class);
        while(df_ex_iterator.hasNext()) {
            Expression expr = df_ex_iterator.next();
            if(expr instanceof ArrayAccess) {
                if(local == 1) {
                    ArrayAccess acc = (ArrayAccess)expr;
                    // if(!(acc.getArrayName() instanceof IDExpression) || (! (primitive.local_buffers.contains((IDExpression)(acc.getArrayName()) ))) )
                    //     continue;

                    if(acc.getArrayName() instanceof IDExpression) {

                        if(primitive.global_buffers.contains((IDExpression)(acc.getArrayName()) )) {
                            if(acc.getNumIndices()>1) {
                                System.out.println("[NOT SUPPORTED] Global array dimension is more than one");
                                continue;
                            }


                            Expression ind_expr = acc.getIndex(0);

                            // System.out.println("Before removing redundants: "+ind_expr);
                            // if(primitive.total_parameter_ind.containsKey((IDExpression)acc.getArrayName())){
                            //     for(Expression check_expr : ind_to_access_info.get(primitive.total_parameter_ind.get((IDExpression)acc.getArrayName())).original_access_map) {
                            //         System.out.println(check_expr);
                            //     }
                            // }
                            

                            int res = removeRedundants(ind_expr, true);

                            // System.out.println("After removing redundants: "+ind_expr);
                            // if(primitive.total_parameter_ind.containsKey((IDExpression)acc.getArrayName())){
                            //     for(Expression check_expr : ind_to_access_info.get(primitive.total_parameter_ind.get((IDExpression)acc.getArrayName())).original_access_map) {
                            //         System.out.println(check_expr);
                            //     }
                            // }
                            if(res==-1){
                                System.out.println("[ERROR] check index access pattern " +ind_expr);
                            }

                            if((ind_expr instanceof IDExpression) && (res==0))
                                ind_expr.swapWith(new IntegerLiteral(0));

                        }
                        else if(primitive.local_buffers.contains((IDExpression)(acc.getArrayName()) )) {
                            acc.addIndex(new NameID("dummy"));
                            Expression ins_ind = new ArrayAccess(new NameID("t_bd"), new NameID("PARTITION_DIM_GRID"));
        //                     Expression ins_ind = new NameID("TASKLET_ID_WG");
                            
                            for(int j=0; j<acc.getNumIndices(); j++) {
                                Expression temp = acc.getIndex(j);
                                acc.setIndex(j, ins_ind);
                                ins_ind = temp;
                            }
                        }
                    }

                    
                }
//                 List<Expression> dim_exprs = acc.getIndices();
//                 dim_exprs.add(0, new NameID("tid"));
                
//                 acc.setIndices(dim_exprs);
            }
        }
        
        /*DepthFirstIterator df_tt_iterator = new DepthFirstIterator(tt);
//         df_ex_iterator.pruneOn(ArrayAccess.class);

        int skip = 0;
        HashSet<Expression> skip_array_access = new HashSet<Expression>();
        
        while(df_tt_iterator.hasNext()) {
//             Expression expr = df_ex_iterator.next();
            Traversable temp_tt = df_tt_iterator.next();
            
            if(!(temp_tt instanceof Expression))
                continue;
            
            Expression expr = (Expression) temp_tt;
            
            if(skip>0) {
                --skip;
                continue;
            }
            if(expr instanceof ArrayAccess) {
                ++skip;
                skip_array_access.add(((ArrayAccess)expr).getArrayName());
                continue;
            }
            
            if(expr instanceof IDExpression) {
                IDExpression id_exp = (IDExpression)expr;
                if(tid_idexp.contains(id_exp) || !(skip_array_access.contains(expr)))
                    id_exp.swapWith(new ArrayAccess(id_exp.clone(), new NameID("tid")));
            }
            
            if((expr instanceof AssignmentExpression) || (expr instanceof BinaryExpression)) {
                
//                 Expression lhs = ((AssignmentExpression)expr).getLHS();
//                 if(lhs instanceof IDExpression) {
//                     if(tid_idexp.contains((IDExpression)lhs))
//                         lhs.swapWith(new ArrayAccess(lhs.clone(), new NameID("tid")));
//                 }
                
                DFIterator<IDExpression> df_id_iter = new DFIterator<IDExpression>(expr, IDExpression.class);
                DFIterator<ArrayAccess> df_arr_iter = new DFIterator<ArrayAccess>(expr, ArrayAccess.class);
                
                HashSet<IDExpression> array_accesses = new HashSet<IDExpression>();
                
                for(ArrayAccess sub_arr :df_arr_iter.getList()) {
                    if(sub_arr.getArrayName() instanceof IDExpression)
                        array_accesses.add((IDExpression)sub_arr.getArrayName());
                }
                
//                 df_id_iter.pruneOn(ArrayAccess.class);
//                 df_id_iter.pruneOn(IDExpression.class);
                
                int skip =0;
                System.out.println("EXPR: " + expr);
                while(df_id_iter.hasNext()) {
                    Expression sub_expr = df_id_iter.next();
//                     if(skip>0) {
//                         if(sub_expr instanceof IDExpression)
//                             --skip; // CAUTION: assuming all array access are using id expressions
//                         else {
//                             System.out.println("ASSUMPTION OF IDExpression followed by ")
//                         }
//                         continue;
//                     }
//                     
//                     if(sub_expr instanceof ArrayAccess) {
//                         skip += 1;
//                     }
                    
                    if(!(sub_expr instanceof IDExpression) || (array_accesses.contains((IDExpression)sub_expr)) || !(tid_idexp.contains((IDExpression)sub_expr)))
                        continue;
                        
                    System.out.println(": " + sub_expr);
                    IDExpression id_exp = (IDExpression)sub_expr;
//                     if(tid_idexp.contains(id_exp))
                    id_exp.swapWith(new ArrayAccess(id_exp.clone(), new NameID("tid")));                    
                    
                }
                    
            }
        }*/
    }
    
    String processProcedureDeclarator(ProcedureDeclarator p_declr) {
        List<Declaration> params_decl = p_declr.getParameters();
        
        String input_0_decl = "#define SIZE_DEFAULT 100\n\n";
                                // "__mram_noinit float GLOBAL_BUFFER[DPU_MRAM_SIZE];\n\n";
        String input_0 = "if(tid == 0) {\n";         //if(INPUT==0) 
        String b_input_0 = "if(TASKLET_ID_WG == 0) {\n";
//         String pre = "DPU_MRAM_HEAP_POINTER";
        String glob_buffer = "DPU_MRAM_HEAP_POINTER";
        String pre = glob_buffer;
        int mram_ptrs = 0;
        
        for(Declaration p_decl: params_decl) {
//             p_declr.removeChild(p_decl);
//             p_decl.detach();
            if(! (p_decl instanceof VariableDeclaration))
                continue;
                
            VariableDeclaration new_p_decl = (VariableDeclaration)p_decl.clone();
            List<Specifier> specs = new_p_decl.getSpecifiers();
//             CodeAnnotation c_annot = new CodeAnnotation("__mram ");
//             c_annot.setPosition(0);
//             new_p_decl.annotate(c_annot);
//             tran_unit.addDeclaration(new_p_decl);
            String append = "";
            String post = "";
            
            boolean pointer_present = false;
            boolean auto_present = false;
            
            int num_declarators = new_p_decl.getNumDeclarators();
            
            String temp = "";
                
            if(new_p_decl.getSpecifiers().contains(Specifier.AUTO)) {
                auto_present = true;
            }
            
            VariableDeclarator v_declr = null;
            
            for(int i=0; i<num_declarators; i++) {
                v_declr = (VariableDeclarator) new_p_decl.getDeclarator(i);
                
                
                if(v_declr.getSpecifiers().contains(PointerSpecifier.UNQUALIFIED)) {
                    pointer_present = true;
                    temp = v_declr.getID().getName();
                    break;
                }
            }
            
//             for(Specifier sp: specs) {
//                 System.out.println(sp.getClass().getName() + " " + sp);
//                 if(sp instanceof PointerSpecifier) {
//                     pointer_present = true;
//                     break;
//                 }
//             }
            String vdecl_str = new_p_decl.toString();
            if(pointer_present) {
                if(!auto_present) {
                    primitive.global_buffers.add(v_declr.getID());
                    
                    append = "__host __dma_aligned __mram_ptr ";
                    // input_0 += temp + " = " + pre;

                    input_0 += temp + " = " + glob_buffer + " + " + temp+"_start; \n";
                    // if(mram_ptrs !=0 ) {
                    //     input_0 += "+"+pre+"_size";    
                    // }
                    // input_0 += ;
                    // input_0 += ";\n";
                    pre = temp;
                    ++mram_ptrs;
                    input_0_decl += "__host __dma_aligned int64_t "+temp + "_start = " + "SIZE_DEFAULT;\n";
                }
                else {
//                     append = "__host __dma_aligned ";
                    vdecl_str = vdecl_str.substring(vdecl_str.indexOf(" ")+1, vdecl_str.length());
                    post = "[N_MULTI_WGS]";
                    b_input_0 += temp + "[t_bd[PARTITION_DIM_GRID]] = " + "mem_alloc(" +temp+ "_size);\n";
                    primitive.local_buffers.add(v_declr.getID());
                    input_0_decl += "__host __dma_aligned int64_t "+temp + "_size = " + "SIZE_DEFAULT;\n";
                }
                
                //input_0_decl += "__host __dma_aligned int "+temp + "_size = " + "SIZE_DEFAULT;\n";
            }
            else {
                primitive.scalars.add(v_declr.getID());
                append = "__host __dma_aligned ";
            }
            
            tran_unit.addDeclarationBefore(proc, createDecls.createCodeAnnotDecl(append + vdecl_str + post + ";"));
            
//             p_decl.setParent(tran_unit);
//             tran_unit.setChild(tran_unit.getChildren().size(), p_decl);
        }
        
        input_0 += "}\n";
        b_input_0 += "}\n"+
                    "barrier_wait(&my_barrier);\n";



        
        tran_unit.addDeclarationBefore(proc, createDecls.createCodeAnnotDecl(input_0_decl));
        return input_0+b_input_0;
    }
    
    boolean swap_statements_IR()
    {
        Iterator<Statement> ref_stmts = ref_swap_stmts.iterator();
        Iterator<Long> swap_dependence_it = swap_dependence.iterator();
        Iterator<Statement> comp_swap_stmts_it = add_comp_swap_stmts.iterator();
        Iterator<Boolean> swap_to_before_it = is_swap_next.iterator();
        Iterator<Integer> skip_swap_stmts_it = skip_swap_stmts.iterator();
        
        Statement prvs_swap_stmt = null;
        long prvs_dependence = 0;
        Statement prvs_ref_stmt = null;
        CompoundStatement prvs_body = null;
        //Iterator<Object> comp_stmts_swap_add = add_comp_swap_stmts.iterator();
        int add_stmts = 0;
        boolean made_not_dep = false;
        int should_skip = 0;
        for(Object st: add_swap_stmts) {
            ++add_stmts;
            int error=0;
            Statement cr_stmt = (Statement)st;
            Statement ref_st = (Statement) ref_stmts.next();
            CompoundStatement c_st = (CompoundStatement) comp_swap_stmts_it.next();
            boolean swap_to_before = swap_to_before_it.next(); 
            int skip = skip_swap_stmts_it.next();
            
            if(should_skip>0) {
                should_skip --;
                continue;
            }
            if(made_not_dep && should_skip==0) {
                swap_to_before = false;
                made_not_dep = false;
            }
                
            long local_dep = swap_dependence_it.next()&partition_idx;
            List<Traversable> pending = new ArrayList<Traversable>(); 
            if(local_dep==0) {
                made_not_dep = true;
                prvs_dependence = 0;
                should_skip += skip;
                
                continue;
            }
//             made_not_dep = false;
            //CompoundStatement comp_st = ((CompoundStatement) comp_stmts_swap_add.next()).clone();
            CompoundStatement body = null;
            
            
            
            if(swap_to_before && (prvs_dependence == local_dep)) {
                if(ref_st == prvs_ref_stmt) {
                    c_st.removeStatement(cr_stmt);
                    //ref_st.swapWith(new_stmt);
                }
                prvs_body.addStatement(cr_stmt);
            }

            else if((IRTools.isAncestorOf(prvs_body, cr_stmt)) && ((prvs_dependence&local_dep)!=0)) {
                body = new CompoundStatement();

                // if(cr_stmt instanceof DeclarationStatement) {
                //     body.addDeclaration(((DeclarationStatement)cr_stmt).getDeclaration().clone());
                // }
                // else
                //     body.addStatement(cr_stmt);

                if(((local_dep^prvs_dependence)&local_dep) == 0)
                    continue;
                
                Statement new_stmt = (Statement)createForLoop(vector_ids, (local_dep^prvs_dependence)&local_dep, (Statement) body);
                
                if(new_stmt==null) {
                    System.out.println("->" + prvs_body +" "+cr_stmt + " "+Long.toBinaryString(local_dep) + " " + Long.toBinaryString(prvs_dependence));
                    System.out.println("[ERROR] NULL");
                }

                ref_st.swapWith(new_stmt);
                body.addStatement(cr_stmt);

            }

            else {
                body = new CompoundStatement();
                // c_st.removeStatement(cr_stmt);

                Statement new_stmt = (Statement)createForLoop(vector_ids, local_dep, (Statement) body);
                
//                 if(swap_to_before) {
//                     c_st.addStatementAfter(prvs_swap_stmt, new_stmt);
//                 }
//                 else {
//                     System.out.println(":: "+ref_st.toString());
                    if(new_stmt==null)
                        System.out.println("[CHECK] NULL");
                    ref_st.swapWith(new_stmt);
                    body.addStatement(cr_stmt);

//                 }
                prvs_body = body;
                prvs_swap_stmt = new_stmt;
            }
            
                //comp_stmt.addStatement((Statement)st);
            
            prvs_dependence = local_dep;
            prvs_ref_stmt = ref_st;
            
        }
        return true;
    }
    
    boolean add_statements_IR() {
        Iterator<Statement> ref_to_stmts_it = ref_to_stmts.iterator();
        Iterator<Long> add_dependence_it = add_dependence.iterator();
        Iterator<Statement> comp_add_stmts_it = add_comp_to_stmts.iterator();
        Iterator<Boolean> add_to_before_it = is_add_next.iterator();
        
        Statement prvs_ref_stmt = null;
        CompoundStatement prvs_body = null;
        long prvs_dependence = 0;
        
        //Iterator<Object> comp_stmts_swap_add = add_comp_swap_stmts.iterator();
        Statement prvs_stmt = null;

        for(Object st: add_to_stmts) {

            int error=0;
            Statement cr_stmt = (Statement)st;
            Statement ref_st = (Statement) ref_to_stmts_it.next();
            CompoundStatement c_st = (CompoundStatement) comp_add_stmts_it.next();
            
            long local_dep = add_dependence_it.next()&partition_idx;
            boolean add_to_before = add_to_before_it.next();
            
            
            if(add_to_before && (prvs_dependence == local_dep) && (ref_st == prvs_ref_stmt)) {

                if(cr_stmt instanceof DeclarationStatement) {
                    prvs_body.addDeclaration(((DeclarationStatement)cr_stmt).getDeclaration().clone());
                }
                else
                    prvs_body.addStatement(cr_stmt);
            }
            else if((IRTools.isAncestorOf(prvs_body, ref_st)) && ((prvs_dependence&local_dep)!=0)) {
                // CAUTION this case not checked
                CompoundStatement body = new CompoundStatement();

                if(cr_stmt instanceof DeclarationStatement) {
                    body.addDeclaration(((DeclarationStatement)cr_stmt).getDeclaration().clone());
                }
                else
                    body.addStatement(cr_stmt.clone());

                body = body.clone();
                Statement new_stmt = (Statement)createForLoop(vector_ids, (local_dep^prvs_dependence)&local_dep, (Statement) body);
                
                if(new_stmt==null) {
                    System.out.println("[CHECK] null");
                    c_st.addStatementAfter(ref_st, cr_stmt);
                }
                else if(add_to_before) {
                    prvs_body.addStatement(new_stmt);
                }
                else {
                    c_st.addStatementAfter(ref_st, new_stmt);
                }

                prvs_body = body;
                prvs_stmt = new_stmt;
            }
            else {
                // if(prvs_body!=null && (IRTools.isAncestorOf(prvs_body, cr_stmt)))
                //     System.out.println("HURRAY"+prvs_body+" " + cr_stmt);
                // if(prvs_body!=null && (IRTools.isAncestorOf(cr_stmt, prvs_body)))
                //     System.out.println("HURRAY2"+prvs_body+" " + cr_stmt);
                // if(cr_stmt instanceof IfStatement) {
                //     System.out.println("KRRR " + cr_stmt);
                // }

                // System.out.println(":( " + cr_stmt);
                // if(prvs_body!=null) {
                //     System.out.println("P BODY: "+ prvs_body);
                // }

                CompoundStatement body = new CompoundStatement();
                //c_st.removeStatement(cr_stmt);
                if(cr_stmt instanceof DeclarationStatement) {
                    body.addDeclaration(((DeclarationStatement)cr_stmt).getDeclaration().clone());
                }
                else
                    body.addStatement(cr_stmt.clone());
                Statement new_stmt = (Statement)createForLoop(vector_ids, local_dep&partition_idx, (Statement) body);
                
                if(add_to_before) {
                    c_st.addStatementAfter(prvs_stmt, new_stmt);
                }
                else {
                    c_st.addStatementAfter(ref_st, new_stmt);

                }
                prvs_body = body;
                prvs_stmt = new_stmt;
                
            }
            
            prvs_dependence = local_dep;
            prvs_ref_stmt = ref_st;

        }
        return true;
    }
    void swap_with_statement(Statement current_stmt, Statement changed_stmt, Statement current_comp_stmt, long dependence_vector, boolean new_one) {

        if(!new_one && ref_swap_stmts.size()!=0)  {
                is_swap_next.add(true);
        }
        else {
            is_swap_next.add(false);
        }
        ref_swap_stmts.add(current_stmt);
        //changed_stmt = if_stmt;
        add_swap_stmts.add(changed_stmt);
        add_comp_swap_stmts.add(current_comp_stmt);
        swap_dependence.add(dependence_vector);
        skip_swap_stmts.add(0);
    }
    
    void add_to_statement(Statement current_stmt, Statement changed_stmt, Statement current_comp_stmt, long dependence_vector, boolean new_one) {

        if(!new_one && ref_to_stmts.size()!=0)  {
                    is_add_next.add(true);

        }
        else {
            is_add_next.add(false);
        }
        //changed_stmt = if_stmt;
        ref_to_stmts.add(current_stmt);
        add_to_stmts.add(changed_stmt);
        add_comp_to_stmts.add(current_comp_stmt);
        add_dependence.add(dependence_vector);
        skip_add_stmts.add(0);
    }
    
    void findArrayAccessPattern(Expression expr) {
        ExpressionStatement expr_stmt = (ExpressionStatement) expr.getStatement();
        System.out.println("expr: " + expr + " |-> " +expr_stmt.getParent().getParent().getClass().getName());
        DFIterator<ArrayAccess> df_iterator = new DFIterator<ArrayAccess>(expr, ArrayAccess.class);

        int is_write = 1;   //read
        if(expr instanceof AssignmentExpression) {
            Expression lhs = ((AssignmentExpression)expr).getLHS();
            if(lhs instanceof ArrayAccess) {
                if(((AssignmentExpression)expr).getOperator()==AssignmentOperator.NORMAL)
                    is_write = 2;   //write
                else
                    is_write = 3;   //both read & write
            }
        }
        
        boolean is_first = true;
        while(df_iterator.hasNext()) {
            ArrayAccess arr_acc = df_iterator.next();
            
            IDExpression id_expr = (IDExpression)arr_acc.getArrayName();

            if(is_first)
                is_first = false;
            else
                is_write = 1;

            if(!param_ind_map.containsKey(id_expr))
                continue;
                // return;
            
            int t_ind = arr_acc.getNumIndices();
            if(t_ind!=1)
                System.out.println("ERROR: Array has mutiple dimensions. Please ensure single indexing");
            Expression ind_expr = arr_acc.getIndex(0);
            
            Expression expr_clone = null;
            /*if(ind_expr instanceof Identifier) {
                Expression var_id = (Identifier)ind_expr;
                Symbol var_sym = ((Identifier)var_id).getSymbol();
//                 if(Definitions.containsKey(var_sym)) {
//                     System.out.println("Expanded: "+ Definitions.get(var_sym).toString());
//                     ((Expression)tt).swapWith(Definitions.get(var_sym).clone());
//                 }
//                 System.out.println("Passed - Expanded: "+ Definitions.get(var_sym).toString());
                expr_clone = Definitions.get(var_sym).clone();
            }
            
            else */
            
            // System.out.println(id_expr.getName()+": "+ ind_expr.toString() );
            // System.out.println(id_expr.getName()+": "+ expr_clone.toString() );
            
            expr_clone = set_coefficients(ind_expr, primitive);
            
            // System.out.println(":: "+id_expr.getName()+": "+ expr_clone.toString() );
                
            if(expr_clone==null) {
                System.out.println("[ERROR]: index analysis not found primitive variables. not consistent " + ind_expr.toString());
                continue;
            }
            if(!param_ind_map.containsKey(id_expr))
                continue;
                
            int p_ind = param_ind_map.get(id_expr);
//             System.out.println(id_expr.toString() + " " + p_ind +" "+ ind_expr.toString() + " " + expr_clone.toString());
            if(ind_to_access_info.get(p_ind)==null) {
                // ind_to_access_map.put(p_ind, new ArrayList<Expression>());
                ind_to_access_info.put(p_ind, new IndexToAccessInfo());
                
                // ind_to_access_type_wmap.put(p_ind, new Long(0));
                // ind_to_access_type_rmap.put(p_ind, new Long(0));
            }
            
//             int nth_time = ind_to_access_map.get(p_ind).size();

            Long added_cond_status = conditionals.getConditionalsToStmt(expr_stmt, primitive, definitions_expand);
            // System.out.println("Conditionals status = " + Long.toBinaryString(added_cond_status));

            IndexToAccessInfo ind_info = ind_to_access_info.get(p_ind); 

            Pair<Integer, Boolean> add_status = ind_info.add_access_expr(ind_expr, expr_clone, is_write, added_cond_status);

            int nth_time = add_status.getFirst();

            // if(add_status.getSecond()) {
            //     Strides ind_coefficients = new Strides();
            //     Expression analyse_expr = Symbolic.simplify(expr_clone.clone());
            //     System.out.println("Analysed Expression = "+analyse_expr);
            //     Pair<Integer, Integer> res_coefficients = AnalyseArrayExpression.analyseIndex(analyse_expr, ind_coefficients, 0, primitive);
            //     if(res_coefficients.getFirst()==-1){
            //         ind_coefficients.is_strided = -1;
            //         System.out.println("[ERROR] check array index for mul variables" + analyse_expr);
            //     }
            //     else {
            //         ind_coefficients.set_offset(analyse_expr);
            //         ind_coefficients.simplify_all();
            //         ind_coefficients.print();

            //         ind_info.set_stride_info(nth_time, ind_coefficients);
            //     }
            // }

            if(nth_time!=0 && (primitive.global_buffers.contains(id_expr))) {
                arr_acc.setIndex(0, new IntegerLiteral(0));
                arr_acc.setIndex(0, new BinaryExpression(ind_expr, BinaryOperator.ADD, new ArrayAccess(new NameID(id_expr.getName()+"_offset"), new IntegerLiteral(nth_time))) );
                // arr_acc.setIndex(0, new BinaryExpression(ind_expr, BinaryOperator.ADD, new NameID(id_expr.getName()+"_offset"+nth_time)) );

                    System.out.println("Before (CHECK 1): ");
                    for(Expression check_expr : ind_info.original_access_map) {
                        System.out.println(check_expr);
                        if(check_expr.equals(ind_expr)) {
                            if(check_expr == ind_expr)
                                System.out.print(" YES ");
                            else
                                System.out.print(" NO ");
                        }
                    }
                    System.out.println(": "+arr_acc.getIndex(0));
            }


            // if(!(ind_to_access_map.get(p_ind).contains(expr_clone))) {
            //     int nth_time = ind_to_access_map.get(p_ind).size();
            //     ind_to_access_map.get(p_ind).add(expr_clone);

            //     // if(is_write) {
            //     //     ind_to_access_type_map.put(p_ind, ind_to_access_type_map.get(p_ind)|(1<<nth_time));
            //     // }
                
            //     if(nth_time!=0 && (primitive.global_buffers.contains(id_expr))) {
            //         arr_acc.setIndex(0, new BinaryExpression(ind_expr.clone(), BinaryOperator.ADD, new NameID(id_expr.getName()+"_offset"+nth_time)) );
            //     }
            // }

            // int nth_time = ind_to_access_map.get(p_ind).indexOf(expr_clone);
            // if(is_write){ 
            //     ind_to_access_type_wmap.put(p_ind, ind_to_access_type_wmap.get(p_ind)|(1<<nth_time));
            // }
            // else {
            //     ind_to_access_type_rmap.put(p_ind, ind_to_access_type_rmap.get(p_ind)|(1<<nth_time));
            // }
//             (ind_to_access_map.get(p_ind)).put(expr_clone, nth_time);
        }
    }
    
	long arrayIndexAnalysis(Traversable obj, Map<Symbol, Expression> pointer_params, Procedure proc) {
		DepthFirstIterator iter = new DepthFirstIterator(obj);
		int count = 0;
		long local_dep = 0;
		while(iter.hasNext()) {
			++count;
			Traversable o = iter.next();
			if(o instanceof ArrayAccess) {
				ArrayAccess acc = (ArrayAccess)o;
				if(!IRTools.containsSymbols(acc, pointer_params.keySet())){
					continue;
				}
				Set<Symbol> acc_symbols = SymbolTools.getAccessedSymbols(acc);
				Symbol array_symbol = SymbolTools.getSymbolOf(acc.getArrayName());

				if(acc.getArrayName() instanceof Identifier && pointer_params.containsKey(array_symbol)) {
					if(acc.getNumIndices() == getEffectiveDimensions(array_symbol)) {
						//System.out.print("index expressions ("+proc.getName().getName()+":"+array_symbol.getSymbolName()+") : ");
						// out.write("index expressions :"+array_symbol.getSymbolName()+") : ");
						List<Expression> ind_exprs = acc.getIndices();
						//long local_dep = 0;
						for(Expression ind_expr: ind_exprs) {
							ind_expr.print(out);
							// out.write(": ");
							Set<Symbol> ind_symbols = SymbolTools.getAccessedSymbols(ind_expr);
							for(Symbol ind_symbol : ind_symbols) {
								//Declaration ind_decl = ind_symbol.getDeclaration();
								//ind_decl.print(out);
								// out.write(ind_symbol.getSymbolName());
								// out.write("  ");
								
								for(int i=0; i<6; i++) {
									if(vector_ids[i]==ind_symbol){
										// out.write(":P ");
										local_dep |= (1<<i);
										//break;
									}
								}
							}	

						}
					}	
													
				}
			}
		}
		// out.write(count + "\n");
		// out.flush();
		
		return local_dep;
	}

	long getDependenceVector(Traversable t, Map<Symbol, Long> dep_vector, Procedure proc) {
		long local_dep = 0;
		Set<Symbol> acc_syms = 	SymbolTools.getAccessedSymbols(t);
		for(Symbol sym : acc_syms) {
			if(dep_vector.containsKey(sym))
                local_dep |= dep_vector.get(sym);

		}
		out.flush();
		return local_dep;
	}

	ForLoop createForLoop(Symbol[] vector_ids, long nesting_idx, Statement for_body) {
		IDExpression iter_id=null;
		
		ForLoop for_loop = null;
		for(long dim = 0; dim<6; dim++  ) {
            long partition_idx = 1<<dim;
            if((nesting_idx&partition_idx) !=0) {
                iter_id = getIteratorVariable(dim);
                Expression limit_id = TranslationVariableNames.getIteratorLimit(dim);
                                                                    
//                 VariableDeclarator declr = new VariableDeclarator(iter_id);
                //Initializer init = new Initializer(new IntegerLiteral(0));
//                 Initializer init = new Initializer(getIteratorStart(dim));
                
//                 declr.setInitializer(init);
//                 VariableDeclaration v_decl = new VariableDeclaration(Specifier.INT, declr);
//                 DeclarationStatement init_iter = new DeclarationStatement(v_decl);

                Expression condition = new BinaryExpression(iter_id.clone(), BinaryOperator.COMPARE_LT, limit_id);
                Expression step = new UnaryExpression(UnaryOperator.POST_INCREMENT, iter_id.clone());
                //Statement body = for_body.clone();
                NullStatement init_iter_null = new NullStatement();        //CHECK
                Statement init_stmt = new ExpressionStatement(new AssignmentExpression(iter_id, AssignmentOperator.NORMAL, TranslationVariableNames.getIteratorStart(dim)));
                for_loop = new ForLoop(init_stmt, condition, step, for_body);
                for_body = for_loop;
            }
        }
		return for_loop;
		
	}
	ArraySpecifier createLimitsArraySpecifier (long access_idx) {    // creates specifier with limits like [ss_tx_end][ss_ty_end]
        List<Expression> indices = new ArrayList<Expression>();
        for(int dim=0; dim<6; dim++) {
            if((access_idx&(1l<<dim))!=0) {
//                 IDExpression acc_id = getIteratorVariable(dim);
                Expression limit_id = new BinaryExpression(TranslationVariableNames.getIteratorLimit(dim), BinaryOperator.SUBTRACT, TranslationVariableNames.getIteratorStart(dim));
                indices.add(limit_id);
            }
//           
        }
        return new ArraySpecifier(indices);
//         v_declrtr.addTrailingSpecifier(new ArraySpecifier(indices));
	}
	
	ArrayAccess createDepArrayAccess(Expression array, long access_idx) {  // creates access with indices like array_name[tx][ty]
        List<Expression> indices = new ArrayList<Expression>();
        for(int dim=0; dim<6; dim++) {
            if((access_idx&(1l<<dim))!=0) {
                Expression acc_id = getIteratorVariable(dim);
                acc_id = new BinaryExpression(acc_id, BinaryOperator.SUBTRACT, TranslationVariableNames.getIteratorStart(dim));
                indices.add(acc_id);
            }
//           
        }
        return new ArrayAccess(array, indices);
	}
	
	long get_ids(Symbol id_sym, IDExpression id_exp, Expression definition) {  //find the work functions and gives the dependence vectors for bare primitive variables
        if(definition instanceof FunctionCall) {
            Expression call_name = ((FunctionCall) definition).getName();
            //out.write("CALL_NAME: "+call_name.toString()+" ");
            if(call_name.toString().equals("get_local_id")) {
                FunctionCall func_call = (FunctionCall) definition;
                int dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
                //out.write("DIM: "+dim+" ");
                vector_ids[dim] = id_sym;
                vector_idexps[dim] = id_exp;

                dep_vector.put(id_sym, 1l<<dim);
                primitive.add_work_item_func(id_exp, dim);
                return (1l<<dim);
            }
            else if(call_name.toString().equals("get_global_id")) {
                FunctionCall func_call = (FunctionCall) definition;
                int dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
                //out.write("DIM: "+dim+" ");
                vector_ids[3+dim] = id_sym;
                vector_idexps[3+dim] = id_exp;

                long g_dep = (1l<<dim)|(1l<<(6+dim));
                dep_vector.put(id_sym, g_dep);
//                 dep_vector.put(id_sym, 1l<<(3+dim));
                // System.out.println("Identified: get_global_id "+dim +" "+ dep_vector.get(id_sym));
                primitive.add_work_item_func(id_exp, 3+dim);
                return g_dep;
            }
            
            else if(call_name.toString().equals("get_group_id")) {
                FunctionCall func_call = (FunctionCall) definition;
                int dim = (int)((IntegerLiteral)func_call.getArgument(0)).getValue();
                //out.write("DIM: "+dim+" ");
                vector_ids[6+dim] = id_sym;
                vector_idexps[6+dim] = id_exp;

                dep_vector.put(id_sym, 1l<<(6+dim));
                // System.out.println("Identified: get_group_id "+dim +" "+ dep_vector.get(id_sym));
                primitive.add_work_item_func(id_exp, 6+dim);
                return (1l<<(6+dim));
            }
            
        }
        return 0l;
	}
	
    boolean getDefinitionsAndDependence(Procedure proc, Traversable trav, boolean new_one) {        // calculates the all variables dependences ang give definitions
        if(trav instanceof DeclarationStatement) {
            Declaration decl = ((DeclarationStatement) trav).getDeclaration();
            boolean set_decl_before = false;
            if(decl instanceof VariableDeclaration) {
                VariableDeclaration v_decl = (VariableDeclaration)decl;
                List<Symbol> decl_syms = v_decl.getDeclaredSymbols();
                List<IDExpression> decl_ids = v_decl.getDeclaredIDs();
                int total = v_decl.getNumDeclarators();
                Iterator<IDExpression> iter_id = decl_ids.iterator();
                Iterator<Symbol> decl_syms_iterator = decl_syms.iterator();
                
                
                for(int index=0; index<total; index++) {
                    Symbol v_sym = decl_syms_iterator.next();
                    VariableDeclarator v_declr = (VariableDeclarator)v_decl.getDeclarator(index);
                    IDExpression v_id = iter_id.next();
                    
//                     if(v_sym==null || v_id==null) {
//                         System.out.println("ERROR: "+v_sym.getSymbolName());
//                     }
//                     else {
//                         System.out.println("NO ERROR: "+v_sym.getSymbolName());
//                     }
                    
                    
                    Initializer init = v_declr.getInitializer();
                    
                    if(init!=null) {
                        //out.write(id_sym.getSymbolName() + " " + init.getClass().getName());
                        Iterator<Traversable> init_value_iterator = init.getChildren().iterator();
                        while(init_value_iterator.hasNext()) {
                            
                            Traversable tt = init_value_iterator.next();
                            //out.write(tt.getClass().getName()+" ");
                            if(tt instanceof Expression) {
                                Definitions.put(v_sym, (Expression)tt);
                                definitions_id.put(v_id, (Expression)tt);
//                                 if( ! primitive.is_iterator_variable(v_id)) {
                                if(!iterator_variables.contains(v_sym)) {
                                    definitions_expand.put(v_id, expand_expression((Expression)tt));
                                }
                                else {
                                    System.out.println("12-05-2021: "+v_id.toString());
                                }
                                long local_dep = getDependenceVector((Expression)tt, dep_vector, proc);
                                dep_vector.put(v_sym, local_dep);
                                
                                boolean loop_tran_check = true;
                                
                                if(tt instanceof FunctionCall) {
                                    local_dep = get_ids(v_sym, v_id, (FunctionCall)tt);
                                    if((local_dep == 0) || (local_dep&(local_dep-1))==0)
                                        loop_tran_check = false;
                                }
                                if(loop_tran_check) {
                                    if((local_dep & partition_idx)!=0) {
                                        VariableDeclarator temp_v_decl = new VariableDeclarator(PointerSpecifier.UNQUALIFIED, v_id.clone());
                                        NestedDeclarator temp_nest_declr = new NestedDeclarator(temp_v_decl);
                                        ArraySpecifier arr_spec = createLimitsArraySpecifier(local_dep&partition_idx);
                                        temp_nest_declr.addTrailingSpecifier(arr_spec);
                                        Expression size_expr = arr_spec.getDimension(0);
                                        for(int ind=1; ind<arr_spec.getNumDimensions(); ind++) {
                                            size_expr = new BinaryExpression(size_expr, BinaryOperator.MULTIPLY, arr_spec.getDimension(ind));
                                        }
                                        size_expr = new BinaryExpression(size_expr, BinaryOperator.MULTIPLY, new SizeofExpression(v_decl.getSpecifiers()));

                                        temp_nest_declr.setInitializer( new Initializer(new FunctionCall(new NameID("mem_alloc"), size_expr)) );
//                                         v_declr.addTrailingSpecifier(createLimitsArraySpecifier(local_dep&partition_idx));
                                        ArrayAccess arr_acc = createDepArrayAccess(new UnaryExpression(UnaryOperator.DEREFERENCE, v_id.clone()), local_dep&partition_idx);
                                        Expression rhs = ((Expression) tt).clone();
                                        ExpressionStatement expr_stmt = createAssignExprStmt(arr_acc, AssignmentOperator.NORMAL, rhs);
//                                         v_declr.setInitializer(null);
                                        v_decl.setChild(index, temp_nest_declr);

                                        // System.out.println("New DECL" + index + ": " + temp_nest_declr);

                                        set_decl_before = true;
                                        //originalExp.add(v_id);
                                        //originalSym.add(v_sym);
                                        originalES.put(v_id, v_sym);
                                        //group_changed_stmts.add(expr_stmt);
                                        //if(group_changed_stmts.size()==1) {
                                        if(tt instanceof FunctionCall) {
//                                                 add_to_statement(cur_stmt, new DeclarationStatement(td_decl), cur_comp_stmt, local_dep, new_one);
//                                             add_to_statement()
                                        }
                                            add_to_statement(cur_stmt, expr_stmt, cur_comp_stmt, local_dep, new_one);
                                            new_one = false;
                                        //}
                                        
                                        
                                    }
                                }
                                break;
                                //return (Expression) tt;
                            }
                            
                            
                        }
                        if(init_value_iterator.hasNext()) {
                            System.out.println("getDefinitions(): May cause inconsistency");
                        }
                    }
                    else {
                        Definitions.put(v_sym, null);
                        definitions_id.put(v_id, null);
                        dep_vector.put(v_sym, 0l);
                    }
                }
            }
            return set_decl_before;
        }
        return false;
    }

	int numSkipStatements (Statement st) {     // number of Statements for a if or for blocks
		int ret_skip = 0;
		DFIterator<Traversable> df_iterator = new DFIterator<Traversable>(st, Statement.class);
// 		if(st instanceof IfStatement) {
            
								//if_df_iterator.pruneOn(Statement.class);
								
//         out.write(df_iterator.next().getClass().getName() + "\n");
//             int temp_skip = 0;
        while(df_iterator.hasNext()) {
            
            Traversable tt = df_iterator.next();
            //out.write(line_num+" "+tt.getClass().getName() + "\n");
            //++line_num;
            if(tt instanceof Statement) {
                //p_iter.next();
                ret_skip++;
            }
            
        }
            
//         }
			/*Statement then_st = ((IfStatement)st).getThenStatement();
			Statement else_st = ((IfStatement)st).getElseStatement();
			if(then_st instanceof CompoundStatement) {
				ret += ((CompoundStatement)then_st).countStatements();			
			}
			if(else_st!=null) {
				if(else_st instanceof CompoundStatement) {
					ret += ((CompoundStatement)else_st).countStatements();			
				}
				//else
					//++ret;
			}*/
// 		}
// 		else if(st instanceof ForLoop) {
//             ret += ((CompoundStatement)((ForLoop)st).getBody()).countStatements();
//         }
//         else if(st instanceof WhileLoop) {
//             ret += ((CompoundStatement)((WhileLoop)st).getBody()).countStatements();
//         }
//         else if(st instanceof DoLoop) {
//             ret += ((CompoundStatement)((DoLoop)st).getBody()).countStatements();
//         }
        out.write("#skipped: "+ret_skip+"\n");
		return ret_skip-1;
			
	}
	
	void createArrayDecl(Symbol array_sym, Expression dimension) {
        ArraySpecifier arr_spec = new ArraySpecifier(dimension);
        IDExpression _id = new Identifier(array_sym);
        VariableDeclarator declr = new VariableDeclarator(_id, arr_spec);
		//Initializer init = new Initializer(new IntegerLiteral(0));
		//declr.setInitializer(init);
		VariableDeclaration v_decl = new VariableDeclaration(Specifier.INT, declr);
		DeclarationStatement init_iter = new DeclarationStatement(v_decl);
        // System.out.println("DECL: "+ v_decl.toString());
	}
	
	Declaration createArrayVariableDeclaration(IDExpression iter_id, List<Expression> init_exprs) {
        VariableDeclarator declr = new VariableDeclarator(iter_id, ArraySpecifier.UNBOUNDED);
        if(init_exprs!=null) {            
            Initializer init = new Initializer(init_exprs);
            declr.setInitializer(init);
        }
        
        VariableDeclaration v_decl = new VariableDeclaration(Specifier.INT, declr);
        return v_decl;
//         DeclarationStatement init_iter = new DeclarationStatement(v_decl);
//         return init_iter;
    }
	
	
	ExpressionStatement createAssignExprStmt(Expression lhs, AssignmentOperator op, Expression rhs) {  //AssignmentOperator.NORMAL
        BinaryExpression b_expr = new AssignmentExpression(lhs, op, rhs);
        return new ExpressionStatement(b_expr);
        
	}
	
	int getEffectiveDimensions(Symbol symbol) {
		int ret = 0;
		if(symbol instanceof NestedDeclarator ||
			symbol.getTypeSpecifiers().contains(PointerSpecifier.UNQUALIFIED))
			ret++;
		List array_specifiers = symbol.getArraySpecifiers();
		if(!array_specifiers.isEmpty() &&
			array_specifiers.get(0) instanceof ArraySpecifier)
			ret += ((ArraySpecifier)array_specifiers.get(0)).getNumDimensions();
		return ret;
	
	}
	
	IDExpression getIteratorVariable(long dim) {
        IDExpression iter_id = null;
        if(dim>9) {
            System.out.println("ERROR: dimension provided is a value not supported");
        }
        if(vector_ids[(int)dim]!=null) {
			iter_id = new Identifier(vector_ids[(int)dim]);
		}
		else {
			String[] id_names = new String[]{"ss_tx", "ss_ty", "ss_tz", "ss_gx", "ss_gy", "ss_gz", "ss_bx", "ss_by", "ss_bz"};
			iter_id = new NameID(id_names[(int)dim]);
		}
           
        return iter_id;
    }
    
    
    Expression set_coefficients(Expression ind_expr, Primitive primitive) {    //CAUTION: pass clone of the expression. //analyzes if all index expressions contains primitves
    
        Expression expr = null;
        if(ind_expr instanceof IDExpression) {
            if(definitions_expand.containsKey((IDExpression)ind_expr)) {
                expr = definitions_expand.get((IDExpression)ind_expr).clone();
            }
            else {
                System.out.println("Definition not found: " + ind_expr);
                expr = ind_expr.clone();
            }
        }
        else {
            expr = expand_expression(ind_expr);
//             expr = ind_expr.clone();
        }
        
//         System.out.println("IND_EXPR: "+ind_expr.toString() + " " + expr.toString());
        
        if(expr instanceof FunctionCall) {
            FunctionCall sub_expr = (FunctionCall)expr;
                    
            int res = primitive.is_primitive(sub_expr);
            
            if(res==-1) {
                return null;
            }
            if(res<64) {
                IDExpression id_iter = TranslationVariableNames.getIteratorVariableHost(res, "", "l");
                
                return id_iter;
            }
        }
        

            DFIterator<Traversable> expr_iterator = new DFIterator<Traversable>(expr);
            expr_iterator.pruneOn(FunctionCall.class);
            expr_iterator.pruneOn(IDExpression.class);
            
            while(expr_iterator.hasNext()) {
                Traversable tt = expr_iterator.next();
                expandBasicExpression(tt, primitive);
            }
//         }
        
        return expr;
    }


    boolean expandBasicExpression(Traversable tt, Primitive primitive)    // if expression primitve then it converts the workfunction calls to host variables
    { 
        if(tt instanceof IDExpression) {
            IDExpression sub_expr = (IDExpression)tt;
            
            int res = primitive.is_primitive(sub_expr);
//                     System.out.println("-- "+ sub_expr.toString()+ " "+res);
            if(res==-1) {
                System.out.println("-1: "+sub_expr.toString());
                return false;
            }
            if(res<64) {
                IDExpression id_iter = TranslationVariableNames.getIteratorVariableHost(res, "", "l");
                
                sub_expr.swapWith(id_iter);
            }
            // else if(res >=96) {
                // sub_expr.swapWith(new NameID("host_"+sub_expr.getName()));
                // sub_expr.swapWith(new NameID(sub_expr.getName()+"l"));
            // }
        }
        else if(tt instanceof FunctionCall) {
            FunctionCall sub_expr = (FunctionCall)tt;
            
            int res = primitive.is_primitive(sub_expr);
            
            if(res==-1) {
                System.out.println("-1: "+sub_expr.toString());
                return false;
            }
            if(res<64) {
                IDExpression id_iter = TranslationVariableNames.getIteratorVariableHost(res, "", "l");
                
                sub_expr.swapWith(id_iter);
//                         System.out.println("-- "+ sub_expr.toString()+ " " +id_iter.toString()+ " "+res);
            }   
        }
        return true;
    }
        Expression expand_expression(Expression expr) {
            List<Traversable> children = expr.getChildren();
            Expression dup_expr = expr.clone();
            
            int child_num = 0;
            for(Traversable child : children) {
                
                if(child instanceof IDExpression) {
                    IDExpression child_id = (IDExpression) child;
                    if(primitive.is_iterator_variable(child_id))  {
//                                 dup_expr.setChild(child_num, new NameID("ss_iter"));
//                             continue;
                    }
                    if(definitions_expand.containsKey(child_id)) {
                        dup_expr.setChild(child_num, definitions_expand.get(child_id).clone());
                    }
                }
                
                else if(!(child instanceof FunctionCall) && !(child instanceof Literal) && (child instanceof Expression)) {
                    
                    Expression sub_expr = (Expression) child;
//                         System.out.println(sub_expr.toString());
                    dup_expr.setChild(child_num, expand_expression(sub_expr));
                }
                ++child_num;
            }
            return dup_expr;
        }
        
        
        Map.Entry<IDExpression, RangeExpression> findLoopIndex(ForLoop L) {
                RangeExpression loop_iter_range = new RangeExpression(new InfExpression(-1), new InfExpression(+1));

                Statement i_stmt = L.getInitialStatement();
                Expression index = null;
                if(i_stmt instanceof ExpressionStatement) {
                    ExpressionStatement init_stmt = (ExpressionStatement)i_stmt;
                    BinaryExpression expr = (BinaryExpression)init_stmt.getExpression();
                    
                    index = expr.getLHS();
                }
                
                else if(i_stmt instanceof DeclarationStatement) {
                    DeclarationStatement init_stmt = (DeclarationStatement)i_stmt;
                    Declaration init_decl = init_stmt.getDeclaration();
                    if(init_decl instanceof VariableDeclaration) {
                        VariableDeclarator var_init_declr = (VariableDeclarator)((VariableDeclaration)init_decl).getDeclarator(0);
                        if(var_init_declr.getInitializer()!=null) {
                            loop_iter_range.setLB((Expression)var_init_declr.getInitializer().getChildren().get(0));
                        }
                    }
                    index = init_decl.getDeclaredIDs().get(0);
                }
                
                if(index == null) {
                    Expression step_expr = L.getStep();
                    System.out.println("STEP EXPR = " + step_expr);
                    if(step_expr instanceof BinaryExpression) {
                        index = ((BinaryExpression)step_expr).getLHS();
                    }
                    else if(step_expr instanceof UnaryExpression) {
                        index = ((UnaryExpression)step_expr).getExpression();
                    }
                }
                
                System.out.println("Iterator variable found: " + index+ " " + i_stmt.getClass().getName());
                
                if(index != null) {
                    Symbol index_sym = SymbolTools.getSymbolOf(index);
                    DFIterator<Statement> df_iterator = new DFIterator<Statement>(L.getBody(), Statement.class);
                    Statement st = null;
                    while(df_iterator.hasNext()) {
                        st = df_iterator.next();
                        if(st instanceof CompoundStatement || st instanceof AnnotationStatement)
                            continue;
                         //rd.getRange(index_sym).toString()
                    }
                    
                    
                    RangeDomain rd = range_map.get(st);
                    
                    if(rd.getRange(index_sym) instanceof RangeExpression) {
                        RangeExpression re= (RangeExpression) rd.getRange(index_sym);
                        if(! (re.getLB() instanceof InfExpression)) {
                            loop_iter_range.setLB(re.getLB().clone());
                        }
                        if(! (re.getUB() instanceof InfExpression)) {
                            loop_iter_range.setUB(re.getUB().clone());
                        }
                    }
                    else if(rd.getRange(index_sym)!=null){
                        loop_iter_range.setUB(rd.getRange(index_sym));
                    }
                    
                    rd = range_map.get(L);
                    if(rd.getRange(index_sym) instanceof RangeExpression) {
                        RangeExpression re= (RangeExpression) rd.getRange(index_sym);
                        if(! (re.getLB() instanceof InfExpression)) {
                            loop_iter_range.setLB(re.getLB().clone());
                        }
                        if(! (re.getUB() instanceof InfExpression)) {
                            loop_iter_range.setUB(re.getUB().clone());
                        }
                    }
                    else if(rd.getRange(index_sym)!=null){
                        loop_iter_range.setLB(rd.getRange(index_sym));
                    }
                    
                    if(loop_iter_range.getLB() instanceof InfExpression) {
                        if(definitions_expand.containsKey((IDExpression)index)) {
                            loop_iter_range.setLB(definitions_expand.get(index));
                        }
                    }
                    
//                     System.out.println(re.getLB().getClass().getName() + " " + re.getLB());
//                     System.out.println(re.getUB().getClass().getName() + " " + re.getUB());

//                     rd = RangeAnalysis.extractRanges(L.getCondition());
                    
                    
                    iterator_variables.add(index_sym);
                    conditionals.addConvertOnes("host_"+((IDExpression)index).getName());
                    definitions_expand.remove(index);
                }
                return new java.util.AbstractMap.SimpleEntry<IDExpression, RangeExpression>((IDExpression)index, loop_iter_range);
        }

	Expression createStepBy2(Expression index) {
                Expression lhs = index.clone();
                Expression rhs = new IntegerLiteral(2);
                Expression indStep2 = new AssignmentExpression(lhs, AssignmentOperator.ADD, rhs);
                return indStep2;
        }

        Expression createIncBy1(Expression index) {
                Expression lhs = index.clone();
                Expression rhs = new IntegerLiteral(1);
                Expression incBy1 = new BinaryExpression(lhs, BinaryOperator.ADD, rhs);
                return incBy1;
        }

        void replaceStep(ForLoop L, Expression new_step) {
                L.setStep(new_step);
        }

        List createNewStmts(ForLoop L) {
                Statement c_stmt = L.getBody();
                List stmts = new ArrayList();

                FlatIterator iter = new FlatIterator(c_stmt);
                while(iter.hasNext()) {
                        Statement st = (Statement)iter.next();
                        stmts.add(st.clone());
                }
                return stmts;
        }

        void insertNewStmts(ForLoop L, List stmts) {
                CompoundStatement c_stmt = (CompoundStatement)L.getBody();

                for(Object st: stmts) {
                        c_stmt.addStatement((Statement)st);
                }
        }

        void replaceExpr(Expression old_expr, Expression new_expr, Traversable t) {
                DepthFirstIterator iter = new DepthFirstIterator(t);
                while(iter.hasNext()) {
                        Object child = iter.next();
                        if(child.equals(old_expr)) {
                                ((Expression)child).swapWith(new_expr.clone());
                        }
                }
        }
}



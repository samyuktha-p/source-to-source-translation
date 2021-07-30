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

public class OpenClDpu extends TransformPass
{
    PrintWriter out;
    public OpenClDpu(Program program) {
            super(program);
            out = new PrintWriter(System.out);
    }

    public String getPassName() {
            return new String("[OpenClDpu]");
    }
    
    // @SuppressWarnings({"all"})
    public void start()
    {
        OneTBtoDPU kernel_transform = new OneTBtoDPU(program);
        HostOneTBtoDPU host_program_transform = new HostOneTBtoDPU(program);
        
        kernel_transform.start();
        
        
        host_program_transform.start();
        
        for(TranslationUnit tu: kernel_transform.tu_kernel_proc) {
            program.addTranslationUnit(tu);
        }

        // Expression example = new BinaryExpression(new NameID("x"), BinaryOperator.ADD, new NameID("y"));
        // Expression example_copy = example;
        // System.out.println("TEST: "+example);
        // DFIterator<IDExpression> df_iter = new DFIterator<IDExpression>(example, IDExpression.class);

        // while(df_iter.hasNext()) {
        //     IDExpression _id = df_iter.next();
        //     _id.swapWith(new NameID("z"));
        // }

        // System.out.println("orig: "+example);
        // System.out.println("copy: "+example_copy);
    
    }

}

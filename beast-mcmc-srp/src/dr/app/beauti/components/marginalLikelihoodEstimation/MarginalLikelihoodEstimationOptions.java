/*
 * MarginalLikelihoodEstimationOptions.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.components.marginalLikelihoodEstimation;

import dr.app.beauti.options.*;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.SequenceErrorType;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MarginalLikelihoodEstimationOptions implements ComponentOptions {

    MarginalLikelihoodEstimationOptions() {
    }

    @Override
    public void createParameters(ModelOptions modelOptions) {
        // nothing to do
    }

    @Override
    public void selectParameters(ModelOptions modelOptions, List<Parameter> params) {
        // nothing to do
    }

    @Override
    public void selectStatistics(ModelOptions modelOptions, List<Parameter> stats) {
        // nothing to do
    }

    @Override
    public void selectOperators(ModelOptions modelOptions, List<Operator> ops) {
        // nothing to do
    }
    //MLE options
    public boolean performMLE = false;
    public int pathSteps = 100;
    public int mleChainLength = 1000000;
    public int mleLogEvery = 1000;
    public String mleFileName = "MLE.log";
    public String pathScheme = "betaquantile";
    public double schemeParameter = 0.30;

}
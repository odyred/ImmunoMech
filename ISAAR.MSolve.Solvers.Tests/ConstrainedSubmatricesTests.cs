﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.LinearSystems;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests
{
    public static class ConstrainedSubmatricesTests
    {
        private static Matrix expectedKff => Matrix.CreateFromArray( new double[,]
        {
            { 9.890109890109889e-01,  0.000000000000000e+00, -3.021978021978022e-01, -1.373626373626375e-02,   1.098901098901099e-01,   0.000000000000000e+00, -2.472527472527472e-01, -1.785714285714285e-01,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  9.890109890109889e-01,  1.373626373626375e-02,  5.494505494505494e-02,   0.000000000000000e+00,  -6.043956043956044e-01, -1.785714285714285e-01, -2.472527472527472e-01,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            {-3.021978021978022e-01,  1.373626373626375e-02,  4.945054945054945e-01, -1.785714285714285e-01,  -2.472527472527472e-01,   1.785714285714285e-01,  5.494505494505494e-02, -1.373626373626375e-02,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            {-1.373626373626375e-02,  5.494505494505494e-02, -1.785714285714285e-01,  4.945054945054945e-01,   1.785714285714285e-01,  -2.472527472527472e-01,  1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 1.098901098901099e-01,  0.000000000000000e+00, -2.472527472527472e-01,  1.785714285714285e-01,   1.978021978021978e+00,   0.000000000000000e+00, -6.043956043956044e-01,  0.000000000000000e+00,  1.098901098901099e-01,   0.000000000000000e+00,  -2.472527472527472e-01, -1.785714285714285e-01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00, -6.043956043956044e-01,  1.785714285714285e-01, -2.472527472527472e-01,   0.000000000000000e+00,   1.978021978021978e+00,  0.000000000000000e+00,  1.098901098901099e-01,  0.000000000000000e+00,  -6.043956043956044e-01,  -1.785714285714285e-01, -2.472527472527472e-01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            {-2.472527472527472e-01, -1.785714285714285e-01,  5.494505494505494e-02,  1.373626373626375e-02,  -6.043956043956044e-01,   0.000000000000000e+00,  9.890109890109889e-01,  0.000000000000000e+00, -2.472527472527472e-01,   1.785714285714285e-01,   5.494505494505494e-02, -1.373626373626375e-02,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            {-1.785714285714285e-01, -2.472527472527472e-01, -1.373626373626375e-02, -3.021978021978022e-01,   0.000000000000000e+00,   1.098901098901099e-01,  0.000000000000000e+00,  9.890109890109889e-01,  1.785714285714285e-01,  -2.472527472527472e-01,   1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   1.098901098901099e-01,   0.000000000000000e+00, -2.472527472527472e-01,  1.785714285714285e-01,  1.978021978021978e+00,   0.000000000000000e+00,  -6.043956043956044e-01,  0.000000000000000e+00,  1.098901098901099e-01,  0.000000000000000e+00, -2.472527472527472e-01, -1.785714285714285e-01,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,  -6.043956043956044e-01,  1.785714285714285e-01, -2.472527472527472e-01,  0.000000000000000e+00,   1.978021978021978e+00,   0.000000000000000e+00,  1.098901098901099e-01,  0.000000000000000e+00, -6.043956043956044e-01, -1.785714285714285e-01, -2.472527472527472e-01,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  -2.472527472527472e-01,  -1.785714285714285e-01,  5.494505494505494e-02,  1.373626373626375e-02, -6.043956043956044e-01,   0.000000000000000e+00,   9.890109890109889e-01,  0.000000000000000e+00, -2.472527472527472e-01,  1.785714285714285e-01,  5.494505494505494e-02, -1.373626373626375e-02,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  -1.785714285714285e-01,  -2.472527472527472e-01, -1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,   1.098901098901099e-01,   0.000000000000000e+00,  9.890109890109889e-01,  1.785714285714285e-01, -2.472527472527472e-01,  1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  1.098901098901099e-01,   0.000000000000000e+00,  -2.472527472527472e-01,  1.785714285714285e-01,  1.978021978021978e+00,  0.000000000000000e+00, -6.043956043956044e-01,  0.000000000000000e+00,  1.098901098901099e-01,   0.000000000000000e+00,  -2.472527472527472e-01, -1.785714285714285e-01 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  -6.043956043956044e-01,   1.785714285714285e-01, -2.472527472527472e-01,  0.000000000000000e+00,  1.978021978021978e+00,  0.000000000000000e+00,  1.098901098901099e-01,  0.000000000000000e+00,  -6.043956043956044e-01,  -1.785714285714285e-01, -2.472527472527472e-01 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00, -2.472527472527472e-01,  -1.785714285714285e-01,   5.494505494505494e-02,  1.373626373626375e-02, -6.043956043956044e-01,  0.000000000000000e+00,  9.890109890109889e-01,  0.000000000000000e+00, -2.472527472527472e-01,   1.785714285714285e-01,   5.494505494505494e-02, -1.373626373626375e-02 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00, -1.785714285714285e-01,  -2.472527472527472e-01,  -1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  1.098901098901099e-01,  0.000000000000000e+00,  9.890109890109889e-01,  1.785714285714285e-01,  -2.472527472527472e-01,   1.373626373626375e-02, -3.021978021978022e-01 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  1.098901098901099e-01,  0.000000000000000e+00, -2.472527472527472e-01,  1.785714285714285e-01,  9.890109890109889e-01,   0.000000000000000e+00,  -3.021978021978022e-01,  1.373626373626375e-02 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00, -6.043956043956044e-01,  1.785714285714285e-01, -2.472527472527472e-01,  0.000000000000000e+00,   9.890109890109889e-01,  -1.373626373626375e-02,  5.494505494505494e-02 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00, -2.472527472527472e-01, -1.785714285714285e-01,  5.494505494505494e-02,  1.373626373626375e-02, -3.021978021978022e-01,  -1.373626373626375e-02,   4.945054945054945e-01,  1.785714285714285e-01 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00, -1.785714285714285e-01, -2.472527472527472e-01, -1.373626373626375e-02, -3.021978021978022e-01,  1.373626373626375e-02,   5.494505494505494e-02,   1.785714285714285e-01,  4.945054945054945e-01 }
        });

        private static Matrix expectedKcf => Matrix.CreateFromArray(new double[,]
        {
            { -3.021978021978022e-01,  -1.373626373626375e-02,  0.000000000000000e+00,  0.000000000000000e+00, -2.472527472527472e-01,  -1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  1.373626373626375e-02,   5.494505494505494e-02,  0.000000000000000e+00,  0.000000000000000e+00, -1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00 },
            { -2.472527472527472e-01,   1.785714285714285e-01,  0.000000000000000e+00,  0.000000000000000e+00, -6.043956043956044e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  -2.472527472527472e-01,  -1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  1.785714285714285e-01,  -2.472527472527472e-01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   1.098901098901099e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00, -2.472527472527472e-01,   1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -6.043956043956044e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  -2.472527472527472e-01,  -1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   1.098901098901099e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  -2.472527472527472e-01,   1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -6.043956043956044e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  -2.472527472527472e-01,  -1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   1.098901098901099e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  -2.472527472527472e-01,   1.785714285714285e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -3.021978021978022e-01,   1.373626373626375e-02,   0.000000000000000e+00,   0.000000000000000e+00 },
            {  0.000000000000000e+00,   0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   1.785714285714285e-01,  -2.472527472527472e-01,   0.000000000000000e+00,   0.000000000000000e+00,  -1.373626373626375e-02,   5.494505494505494e-02,   0.000000000000000e+00,   0.000000000000000e+00 }
        });

        private static Matrix expectedKcc => Matrix.CreateFromArray(new double[,]
        {
            { 4.945054945054945e-01,  1.785714285714285e-01,  5.494505494505494e-02,  1.373626373626375e-02,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 },
            { 1.785714285714285e-01,  4.945054945054945e-01, -1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 },
            { 5.494505494505494e-02, -1.373626373626375e-02,  9.890109890109889e-01,  0.000000000000000e+00,  5.494505494505494e-02,  1.373626373626375e-02,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 },
            { 1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  9.890109890109889e-01, -1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  5.494505494505494e-02, -1.373626373626375e-02,  9.890109890109889e-01,  0.000000000000000e+00,  5.494505494505494e-02,  1.373626373626375e-02,  0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  9.890109890109889e-01, -1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  0.000000000000000e+00 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  5.494505494505494e-02, -1.373626373626375e-02,  9.890109890109889e-01,  0.000000000000000e+00,  5.494505494505494e-02,  1.373626373626375e-02 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  1.373626373626375e-02, -3.021978021978022e-01,  0.000000000000000e+00,  9.890109890109889e-01, -1.373626373626375e-02, -3.021978021978022e-01 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  5.494505494505494e-02, -1.373626373626375e-02,  4.945054945054945e-01, -1.785714285714285e-01 },
            { 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,  1.373626373626375e-02, -3.021978021978022e-01, -1.785714285714285e-01,  4.945054945054945e-01 }
        });

        private static Matrix expectedCondensedK => Matrix.CreateFromArray(new double[,]
        {
            {  1.899630077471670e-01,   8.935750275975297e-02,  -2.292666588239909e-01,   2.466969780372398e-02, -5.106070000547310e-02,  -4.103313066810881e-02,   3.006934549425098e-02,  -4.314166600524617e-02,   6.029500558804593e-02,  -2.985240389012186e-02 },
            {  8.935750275975296e-02,   4.195638480453572e-01,  -1.080162197198449e-01,  -3.359247222029917e-01, -2.218766494944725e-02,  -4.255213035151419e-02,   1.099397801941727e-02,  -2.610167652330916e-02,   2.985240389012187e-02,  -1.498531896754213e-02 },
            { -2.292666588239909e-01,  -1.080162197198449e-01,   4.704867695897013e-01,   3.179465334271065e-02, -2.231042172131791e-01,   8.206626133621771e-02,  -4.818523904678258e-02,   5.149283060333881e-03,   3.006934549425101e-02,  -1.099397801941726e-02 },
            {  2.466969780372398e-02,  -3.359247222029917e-01,   3.179465334271064e-02,   7.777399531778657e-01, -9.445673409134700e-02,  -3.641437143084315e-01,  -5.149283060333874e-03,  -5.156984014313339e-02,   4.314166600524619e-02,  -2.610167652330915e-02 },
            { -5.106070000547312e-02,  -2.218766494944726e-02,  -2.231042172131790e-01,  -9.445673409134700e-02,  5.483298344373040e-01,   0.000000000000000e+00,  -2.231042172131791e-01,   9.445673409134701e-02,  -5.106070000547303e-02,   2.218766494944724e-02 },
            { -4.103313066810881e-02,  -4.255213035151421e-02,   8.206626133621772e-02,  -3.641437143084315e-01, -1.040834085586084e-17,   8.133916893198914e-01,  -8.206626133621771e-02,  -3.641437143084315e-01,   4.103313066810880e-02,  -4.255213035151419e-02 },
            {  3.006934549425100e-02,   1.099397801941727e-02,  -4.818523904678257e-02,  -5.149283060333854e-03, -2.231042172131791e-01,  -8.206626133621769e-02,   4.704867695897013e-01,  -3.179465334271068e-02,  -2.292666588239909e-01,   1.080162197198449e-01 },
            { -4.314166600524617e-02,  -2.610167652330916e-02,   5.149283060333874e-03,  -5.156984014313338e-02,  9.445673409134701e-02,  -3.641437143084315e-01,  -3.179465334271069e-02,   7.777399531778657e-01,  -2.466969780372403e-02,  -3.359247222029917e-01 },
            {  6.029500558804596e-02,   2.985240389012189e-02,   3.006934549425102e-02,   4.314166600524620e-02, -5.106070000547307e-02,   4.103313066810883e-02,  -2.292666588239910e-01,  -2.466969780372398e-02,   1.899630077471670e-01,  -8.935750275975297e-02 },
            { -2.985240389012186e-02,  -1.498531896754213e-02,  -1.099397801941726e-02,  -2.610167652330916e-02,  2.218766494944726e-02,  -4.255213035151419e-02,   1.080162197198449e-01,  -3.359247222029917e-01,  -8.935750275975302e-02,   4.195638480453572e-01 }
        });

        [Fact]
        public static void CheckStaticCondensation()
        {
            double tolerance = 1E-10;

            Model model = CreateModel();
            int id = model.Subdomains.First().ID;

            var solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var problem = new ProblemStructural(model, solver);

            // Prepare model
            model.ConnectDataStructures();

            // Order dofs and initialize linear system
            solver.OrderDofs(true);
            ILinearSystem linearSystem = solver.LinearSystems.First().Value;
            linearSystem.Reset(); // Necessary to define the linear system's size 

            // Build and assign global matrices
            (IMatrixView Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) = 
                problem.CalculateSubMatrices(model.Subdomains.First());
            linearSystem.Matrix = Kff;

            // Static condensation: Kcondensed = Kcc - Kcf * inv(Kff) * Kfc
            Dictionary<int, Matrix> invKffTimesKfc = solver.InverseSystemMatrixTimesOtherMatrix(
                new Dictionary<int, IMatrixView>() { { id, Kfc } });
            IMatrixView condensedK = Kcc.Subtract(Kcf.MultiplyRight(invKffTimesKfc[id]));

            // Checks
            Assert.True(expectedCondensedK.Equals(condensedK, tolerance));
        }

        [Fact]
        public static void CheckSubmatrices()
        {
            double tolerance = 1E-10;

            Model model = CreateModel();
            int id = model.Subdomains.First().ID;

            var solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var problem = new ProblemStructural(model, solver);

            // Prepare model
            model.ConnectDataStructures();

            // Order dofs and initialize linear system
            solver.OrderDofs(true);
            ILinearSystem linearSystem = solver.LinearSystems.First().Value;
            linearSystem.Reset(); // Necessary to define the linear system's size 

            // Build and assign global matrices
            (IMatrixView Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) = 
                problem.CalculateSubMatrices(model.Subdomains.First());

            // Checks
            Assert.True(expectedKff.Equals(Kff, tolerance));
            Assert.True(expectedKcf.Transpose().Equals(Kfc, tolerance));
            Assert.True(expectedKcf.Equals(Kcf, tolerance));
            Assert.True(expectedKcc.Equals(Kcc, tolerance));

        }

        private static Model CreateModel()
        {
            int numElementsX = 2, numElementsY = 4; // for the hardcoded stiffness matrices

            var cantileverBuilder = new CantileverBeam.Builder();
            cantileverBuilder.Length = numElementsX;
            cantileverBuilder.Height = numElementsY;
            cantileverBuilder.Width = 1.0;
            cantileverBuilder.YoungModulus = 1.0;
            cantileverBuilder.PoissonRatio = 0.3;
            CantileverBeam cantilever = cantileverBuilder.BuildWithQuad4Elements(numElementsX, numElementsY);

            return cantilever.Model;
        }
    }
}

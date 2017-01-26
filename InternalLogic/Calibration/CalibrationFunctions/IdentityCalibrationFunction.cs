﻿using System.Text;

namespace InternalLogicCalibration
{
    public class IdentityCalibrationFunction : CalibrationFunction
    {

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Identity");
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            return 0;
        }

        #endregion Internal Methods

    }
}
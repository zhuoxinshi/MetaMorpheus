﻿using EngineLayer;
using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using Spectra;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {

        #region Private Fields

        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForCalibrationTask> ModFileListInWindow = new ObservableCollection<ModListForCalibrationTask>();

        #endregion Private Fields

        #region Public Constructors

        public CalibrateTaskWindow(ObservableCollection<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices(modList);

            TheTask = new CalibrationTask(modList);
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Calibration Task";
        }

        public CalibrateTaskWindow(CalibrationTask myCalibrateTask, ObservableCollection<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices(modList);

            TheTask = myCalibrateTask;
            UpdateFieldsFromTask(TheTask);
        }

        #endregion Public Constructors

        #region Internal Properties

        internal CalibrationTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        private void UpdateFieldsFromTask(CalibrationTask task)
        {
            missedCleavagesTextBox.Text = task.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.Protease;
            maxModificationIsoformsTextBox.Text = task.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.InitiatorMethionineBehavior;

            productMassToleranceTextBox.Text = task.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = (int)task.ProductMassTolerance.Unit;
            precursorMassToleranceTextBox.Text = task.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = (int)task.PrecursorMassTolerance.Unit;

            bCheckBox.IsChecked = task.BIons;
            yCheckBox.IsChecked = task.YIons;
            for (int i = 0; i < ModFileListInWindow.Count; i++)
            {
                if (task.ListOfModListsForCalibration[i].Fixed)
                    ModFileListInWindow[i].Fixed = true;
                if (task.ListOfModListsForCalibration[i].Variable)
                    ModFileListInWindow[i].Variable = true;
                if (task.ListOfModListsForCalibration[i].Localize)
                    ModFileListInWindow[i].Localize = true;
            }
            modificationsDataGrid.Items.Refresh();
        }

        private void PopulateChoices(ObservableCollection<ModList> modList)
        {
            foreach (Protease protease in ProteaseDictionary.Instance.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                productMassToleranceComboBox.Items.Add(toleranceUnit);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                precursorMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            foreach (var uu in modList)
                ModFileListInWindow.Add(new ModListForCalibrationTask(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;
        }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.Protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;

            TheTask.BIons = bCheckBox.IsChecked.Value;
            TheTask.YIons = yCheckBox.IsChecked.Value;

            TheTask.ListOfModListsForCalibration = ModFileListInWindow.ToList();


            TheTask.ProductMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.ProductMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;


            TheTask.PrecursorMassTolerance.Value = double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.PrecursorMassTolerance.Unit = (ToleranceUnit)precursorMassToleranceComboBox.SelectedIndex;


            DialogResult = true;
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            var hm = ye.Content as TextBlock;
            if (hm != null && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        #endregion Private Methods

    }
}
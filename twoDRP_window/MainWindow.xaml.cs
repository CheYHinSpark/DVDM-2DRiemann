using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace twoDRP_window
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            string compute = computeSetTxt.Text;
            string init = initSetTxt.Text;
            string mode = modeBox.SelectedIndex.ToString();
            try
            {
                Process.Start(@".\twoD_Riemann.exe", "c:" + compute + " i:" + init + " m:" + mode);
            }
            catch
            {
                MessageBox.Show("启动失败！");
            }
        }
    }
}

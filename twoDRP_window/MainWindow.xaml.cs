using System.Diagnostics;
using System.Windows;

namespace twoDRP_window
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
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

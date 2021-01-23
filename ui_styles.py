########################################################################################################################
# BY: HALIL IBRAHIM OZDEMIR
# PROJECT MADE WITH: Qt Designer and PyQt5
# V: 0.1
########################################################################################################################

class Style:

    style_bt_standard = (
    """
    QPushButton {
        font-size: 10pt; font-family: Segoe UI;
        font: bold;
        background-image: ICON_REPLACE;
        background-position: left center;
        background-repeat: no-repeat;
        border: none;
        border-left: 28px solid rgb(27, 29, 35);
        background-color: rgb(27, 29, 35);
        text-align: left;
        padding-left: 45px;
        margin-left:-3px;
    }
    QPushButton[Active=true] {
        background-image: ICON_REPLACE;
        background-position: left center;
        background-repeat: no-repeat;
        border: none;
        border-left: 28px solid rgb(27, 29, 35);
        border-right: 5px solid rgb(44, 49, 60);
        background-color: rgb(27, 29, 35);
        text-align: left;
        padding-left: 45px;
    }
    QPushButton:hover {
        background-color: rgb(33, 37, 43);
        border-left: 28px solid rgb(33, 37, 43);
    }
    QPushButton:pressed {
        background-color: rgb(85, 170, 255);
        border-left: 28px solid rgb(85, 170, 255);
    }
    
    
    """
    )

    MessageBox_stylesheet = """

            QWidget
            {
                color: white;
                background-color: #323232;
                selection-background-color:#323232;
                selection-color: black;
                background-clip: border;
                border-image: none;
                border: 0px transparent black;
                outline: 0;
                
                border-radius: 12px;
                border: 1.5px solid;
                padding: 6px;
                border-color: rgb(255, 170, 0);

            }
            
            QLabel
            {
                font: bold 12px;
                border: none;
            }

            QPushButton 
            {
                color: white; 
                font-weight: bold; 
                font-size: 10px; 
                border: 2px solid rgb(52, 59, 72); 
                border-radius: 5px; 
                background-color:  rgb(22, 200, 244); 
                margin-top:1px; 
                margin-bottom: 1px; 
                border-width: 1px; 
                padding: 5px; 
                outline: none;
            }

            QPushButton:hover 
            { 
                background-color: rgb(255, 17, 100); 
                border: 2px solid rgb(61, 70, 86);
            }

            QPushButton:pressed 
            { 
                background-color:  rgb(15, 133, 163); 
                border: 2px solid rgb(43, 50, 61);
            }


        """

    MessageBox_critical_stylesheet = """

                QWidget
                {
                    color: white;
                    background-color: #323232;
                    selection-background-color:#323232;
                    selection-color: black;
                    background-clip: border;
                    border-image: none;
                    border: 0px transparent black;
                    outline: 0;

                    border-radius: 12px;
                    border: 1.5px solid;
                    padding: 6px;

                }

                QLabel
                {
                    font: bold 12px;
                    border: none;
                }

            """

    QToolTip_stylesheet = """
    
                QToolTip 
                { 
                    background-color: black; 
                    color: white; 
                    background-clip: border;
                    border-image: none;
                    border: 0px transparent black;
                    outline: 0;
                    border: black solid 1px;
                    padding: 4px;
                    border-radius: 6px;
                    border: 1px solid;
                    font: bold 10px;
                    border-color:  rgb(157, 90, 198);

                }
    
            """


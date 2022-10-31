import smtplib
from email.utils import formataddr
from email.message import EmailMessage
from email.mime.application import MIMEApplication
from email.mime.text import MIMEText
from os.path import basename
from email.mime.multipart import MIMEMultipart


def sendMailOutlook(mss="Tu script ha finalizado de correr", secrests_path='secrets.txt', files=None, fail=False):
    """Envia un correo electrónico a travez de Outlook/Hotmial/Live. Se debe crear un archivos secrets.txt que contenga:

    Primera linea: correo electronico del que envia el mensaje. Ejemplo: espaiderman@hotmail.com
    Segunda linea: contraseña del correo electrónico
    Tercera linea: correo electrónico del destinatario. Ejemplo: fatman@gmail.com

    El archivo completo se ve así:

        espaiderman@hotmail.com
        MeDanMiedoLosGatos
        fatman@gmail.com

    Args:
        mss (str, optional): Cuerpo del correo electrónico. Defaults to "Tu script ha finalizado de correr".
        secrests_path (str, optional): Ruta al archivo secrets.txt . Defaults to 'secrets.txt'.
    """

    with open(secrests_path) as f:
        e, p, r = f.readlines()
        e = e.replace("\n", "")
        p = p.replace("\n", "")
        r = r.replace("\n", "")
        y = {"sender_email": e, "password": p, "recipient_email": r}

    msg = MIMEMultipart()
    msg['From'] = formataddr(('Python' + fail*'FAIL!!!', y["sender_email"]))
    msg['To'] = y['recipient_email']
    msg['Subject'] = "Tu script ha" + fail * \
        ' fallado!!!' + (1-fail)*' finalizado!'
    msg.attach(MIMEText(mss))
    # 1/0
    for f in files or []:
        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=basename(f)
            )
        # After the file is closed
        part['Content-Disposition'] = 'attachment; filename="%s"' % basename(f)
        msg.attach(part)
    s = smtplib.SMTP("smtp-mail.outlook.com", 587)
    # Hostname to send for this command defaults to the fully qualified domain name of the local host.
    s.ehlo()
    s.starttls()  # Puts connection to SMTP server in TLS mode
    s.ehlo()
    s.login(y["sender_email"], y['password'])
    s.sendmail(y["sender_email"], y["recipient_email"], msg.as_string())
    s.quit()


def sendMailGmail(mss="Tu script ha finalizado de correr", secrests_path='secrets.txt'):
    """Envia un correo electrónico a travez de Gmail. Se debe crear un archivos secrets.txt que contenga:

    Primera linea: correo electronico del que envia el mensaje. Ejemplo: espaiderman@gmail.com
    Segunda linea: contraseña del correo electrónico
    Tercera linea: correo electrónico del destinatario. Ejemplo: fatman@gmail.com

    El archivo completo se ve así:

        espaiderman@gmail.com
        MeDanMiedoLosGatos
        fatman@gmail.com

    Args:
        mss (str, optional): Cuerpo del correo electrónico. Defaults to "Tu script ha finalizado de correr".
        secrests_path (str, optional): Ruta al archivo secrets.txt . Defaults to 'secrets.txt'.
    """
    with open(secrests_path) as f:
        e, p, r = f.readlines()
        e = e.replace("\n", "")
        p = p.replace("\n", "")
        r = r.replace("\n", "")
        y = {"sender_email": e, "password": p, "recipient_email": r}

    msg = EmailMessage()
    msg['From'] = formataddr(('Python', y["sender_email"]))
    msg['To'] = y['recipient_email']
    msg['Subject'] = "Tu script ha finalizado!"
    msg.set_content(mss)
    # 1/0
    s = smtplib.SMTP("smtp.gmail.com", 587)
    # Hostname to send for this command defaults to the fully qualified domain name of the local host.
    s.ehlo()
    s.starttls()  # Puts connection to SMTP server in TLS mode
    s.ehlo()
    s.login(y["sender_email"], y['password'])
    s.sendmail(y["sender_email"], y["recipient_email"], msg.as_string())
    s.quit()


if __name__ == '__main__':
    sendMailOutlook()

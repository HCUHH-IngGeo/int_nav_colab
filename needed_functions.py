import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pyproj import Proj, transform
from sympy import re, sqrt
from colour import Color
import numpy as np
from datetime import timedelta

mapbox_access_token = "pk.eyJ1IjoianN0ODkiLCJhIjoiY2sybjN0b2w2MG1tYjNjcGV5eGtlYjg5NyJ9.A9sSZPYzTYJL1YErdmwoKA"


def load_data(file, skiprows, transform_to_utm=True):
    # read the csv file and skip rows correspondig to inpunt
    df = pd.read_csv(file, skiprows=skiprows)
    # renaming some colums to get easier access on columns
    df.rename(columns={df.columns[0]: "GPST", df.columns[1]: "lat", df.columns[2]: "lon", df.columns[3]: 'alt',
                       df.columns[4]: 'Q', df.columns[5]: 'ns', df.columns[6]: 'sdn', df.columns[7]: 'sde',
                       df.columns[8]: 'sdu', df.columns[9]: 'sdne', df.columns[10]: 'sdeu', df.columns[11]: 'sdun',
                       df.columns[12]: 'age', df.columns[13]: 'ratio'}, inplace=True)
    # change datatype from string to timestamp by given format
    df.GPST = pd.to_datetime(df.GPST, format="%Y/%m/%d %H:%M:%S.%f")
    df.sde *= 3
    df.sdn *= 3


    # since the coordinates are located in wgs85 (lat,lon,alt) we want to transform them to utm
    if transform_to_utm:
        inProj = Proj('epsg:4326')  # projection for wgs84
        outProj = Proj('epsg:25832')  # projection for utm zone 32
        l, l2 = transform(inProj, outProj, list(df.lat),
                          list(df.lon))  # use the transform function to transform the coordinates as lists
        df['x'] = l  # replace the longitude with the eastern values
        df['y'] = l2  # replace the latitude with the northern values
    # print(df.head())
    df.x = df.x.round(3)
    df.y = df.y.round(3)

    df = add_uncertainty(df)
    df = transform_uncertainty(df)

    df = df.set_index(df.GPST)
    df = df.resample('1000ms').mean().interpolate()
    df['GPST']=df.index
    df=df.reset_index(drop=True)

    return df


def get_quality_level(df):
    df_q1 = len(df[df.Q == 1]) / len(df) * 100
    df_q2 = len(df[df.Q == 2]) / len(df) * 100
    df_q5 = len(df[df.Q == 5]) / len(df) * 100
    return df_q1, df_q2, df_q5


# to get the intersectionpoints easternvalue, we have to insert the curve of a first order (y=m*x+b) into the normal
# ellipseequation ((x−x0)**2/a**2+(y−y0)**2/b**2=1) and solve the euqation to x, therefor we get two solutions
def getX(a, b, x1, x2, y1, y2):
    rw1 = -a * b * (y1 - y2) / sqrt(
        a ** 2 * x1 ** 2 - 2 * a ** 2 * x1 * x2 + a ** 2 * x2 ** 2 + b ** 2 * y1 ** 2 - 2 * b ** 2 * y1 * y2 + b ** 2 * y2 ** 2) + x1
    rw2 = a * b * (y1 - y2) / sqrt(
        a ** 2 * x1 ** 2 - 2 * a ** 2 * x1 * x2 + a ** 2 * x2 ** 2 + b ** 2 * y1 ** 2 - 2 * b ** 2 * y1 * y2 + b ** 2 * y2 ** 2) + x1
    return float(re(rw1)), float(re(rw2))


# same as above, but taken the eastern values into the normal ellipse qiation, to get the northern values
def getY(a, b, x1, y1, rw1, rw2):
    hw1 = (a * y1 - b * sqrt((a - rw1 + x1) * (a + rw1 - x1))) / a
    hw2 = (a * y1 + b * sqrt((a - rw2 + x1) * (a + rw2 - x1))) / a
    return float(re(hw1)), float(re(hw2))


def add_uncertainty(df):
    # since we know the (aproximatily) direction of the trajecory, we want to draw the uncertaity besides the trajecotry
    # to get these, we can make use of the std in eastern and northern (which are given by rtklib) we can calculate the intersection
    # of the resulting ellipse with the perpendicular line to the trajectory in the first point of a trajectory segement.

    l1 = []
    l2 = []
    for index, row in df.iterrows():
        if index == len(df) - 1:
            x1 = row.x
            y1 = row.y
            x2 = df.x.iloc[0]
            y2 = df.y.iloc[0]
        else:
            x1 = row.x
            y1 = row.y
            x2 = df.x.iloc[index + 1]
            y2 = df.y.iloc[index + 1]

        a = row.sde
        b = row.sdn

        rw1, rw2 = getX(a, b, x1, x2, y1, y2)
        hw1, hw2 = getY(a, b, x1, y1, rw1, rw2)

        l1.append([float(rw1), float(hw1)])
        l2.append([float(rw2), float(hw2)])

    df_l1 = pd.DataFrame(l1, columns=['x', 'y'])
    df_l1 = df_l1.round(3)
    df_l2 = pd.DataFrame(l2, columns=['x', 'y'])
    df_l2 = df_l2.round(3)

    df['l_x'] = df_l1.x
    df['l_y'] = df_l1.y
    df['r_x'] = df_l2.x
    df['r_y'] = df_l2.y

    return df
    # figure.add_trace(go.Scatter(x=df_l1.x, y=df_l1.y, line=dict(color=color, dash='dash')))
    # figure.add_trace(go.Scatter(x=df.x, y=df.y, line=dict(color=color)))
    # figure.add_trace(go.Scatter(x=df_l2.x, y=df_l2.y, line=dict(color=color, dash='dash')))


def transform_uncertainty(df):
    inProj = Proj('epsg:25832')  # projection for wgs84
    outProj = Proj('epsg:4326')  # projection for utm zone 32
    l, l2 = transform(inProj, outProj, list(df.l_x),
                      list(df.l_y))  # use the transform function to transform the coordinates as lists
    df['l_lat'] = l  # replace the latit
    df['l_lon'] = l2  # replace the longitude with the eastern values

    r, r2 = transform(inProj, outProj, list(df.r_x),
                      list(df.r_y))  # use the transform function to transform the coordinates as lists

    df['r_lat'] = r  # replace the latit
    df['r_lon'] = r2  # replace the longitude with the eastern values
    return df


def add_trace_and_uncertainty(df, figure, color_of_trajectory, name):
    uncertainty_color = Color(color_of_trajectory)
    uncertainty_color.set_hue(uncertainty_color.get_hue() + .15)

    figure.add_trace(go.Scattermapbox(
        lat=df.lat,
        lon=df.lon,
        mode='markers+lines',
        name=f'{name}',
        line={'color': color_of_trajectory}
    ))

    figure.add_trace(go.Scattermapbox(
        lat=df.l_lat,
        lon=df.l_lon,
        mode='lines',
        name=f'{name} - left uncertainty',
        line={'color': f'{uncertainty_color.hex}'}
    ))

    figure.add_trace(go.Scattermapbox(
        lat=df.r_lat,
        lon=df.r_lon,
        mode='lines',
        name=f'{name} - right uncertainty',
        line={'color': f'{uncertainty_color.hex}'}
    ))


def xcorr(x_red, y_red, normed=True, maxlags=None):
    Nx = len(x_red)
    if Nx != len(y_red):
        raise ValueError('x and y must be equal length')

    c = np.correlate(x_red, y_red, mode='full')

    if normed:
        c /= np.sqrt(np.dot(x_red, x_red) * np.dot(y_red, y_red))

    if maxlags is None:
        maxlags = Nx - 1

    if maxlags >= Nx or maxlags < 1:
        raise ValueError('maglags must be None or strictly '
                         'positive < %d' % Nx)

    lags = np.arange(-maxlags, maxlags + 1)
    c = c[Nx - 1 - maxlags:Nx + maxlags]

    return lags, c  # , a, b


def plotData(data):
    figure = make_subplots(2, 1, shared_xaxes=True)
    figure.add_trace(go.Scatter(
        x=data.index,
        y=data.ACC_X,
        mode='lines',
        name=f'Acceleration X',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data.index,
        y=data.ACC_Y,
        mode='lines',
        name=f'Acceleration Y',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data.index,
        y=data.ACC_Z,
        mode='lines',
        name=f'Acceleration Z',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data.index,
        y=data.GYR_X,
        mode='lines',
        name=f'Rotation X',
    ), row=2, col=1)

    figure.add_trace(go.Scatter(
        x=data.index,
        y=data.GYR_Y,
        mode='lines',
        name=f'Rotation Y',
    ), row=2, col=1)

    figure.add_trace(go.Scatter(
        x=data.index,
        y=data.GYR_Z,
        mode='lines',
        name=f'Rotation Z',
    ), row=2, col=1)
    return figure


def plotAccData(data):
    figure = make_subplots(3, 1, shared_xaxes=True)
    figure.add_trace(go.Scatter(
        x=data.TIME,
        y=data.ACC_X,
        mode='lines',
        name=f'Acceleration X',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data.TIME,
        y=data.ACC_Y,
        mode='lines',
        name=f'Acceleration Y',
    ), row=2, col=1)

    figure.add_trace(go.Scatter(
        x=data.TIME,
        y=data.ACC_Z,
        mode='lines',
        name=f'Acceleration Z',
    ), row=3, col=1)
    return figure


def plotAcc_Kalib(data1, data2, data3, data4, data5, data6):
    figure = make_subplots(3, 2)

    figure.add_trace(go.Scatter(
        x=data1.TIME,
        y=data1.ACC_X,
        mode='lines',
        name=f'ACC X where X = 9.81',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data1.TIME,
        y=data1.ACC_Y,
        mode='lines',
        name=f'ACC Y where X = 9.81',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data1.TIME,
        y=data1.ACC_Z,
        mode='lines',
        name=f'ACC Z where X = 9.81',
    ), row=1, col=1)

    figure.add_trace(go.Scatter(
        x=data2.TIME,
        y=data2.ACC_X,
        mode='lines',
        name=f'ACC X where X = -9.81',
    ), row=1, col=2)

    figure.add_trace(go.Scatter(
        x=data2.TIME,
        y=data2.ACC_Y,
        mode='lines',
        name=f'ACC Y where X = -9.81',
    ), row=1, col=2)

    figure.add_trace(go.Scatter(
        x=data2.TIME,
        y=data2.ACC_Z,
        mode='lines',
        name=f'ACC Z where X = -9.81',
    ), row=1, col=2)

    figure.add_trace(go.Scatter(
        x=data3.TIME,
        y=data3.ACC_X,
        mode='lines',
        name=f'ACC X where Y = 9.81',
    ), row=2, col=1)

    figure.add_trace(go.Scatter(
        x=data3.TIME,
        y=data3.ACC_Y,
        mode='lines',
        name=f'ACC Y where Y = 9.81',
    ), row=2, col=1)

    figure.add_trace(go.Scatter(
        x=data3.TIME,
        y=data3.ACC_Z,
        mode='lines',
        name=f'ACC Z where Y = 9.81',
    ), row=2, col=1)

    figure.add_trace(go.Scatter(
        x=data4.TIME,
        y=data4.ACC_X,
        mode='lines',
        name=f'ACC X where Y = -9.81',
    ), row=2, col=2)

    figure.add_trace(go.Scatter(
        x=data4.TIME,
        y=data4.ACC_Y,
        mode='lines',
        name=f'ACC Y where Y = -9.81',
    ), row=2, col=2)

    figure.add_trace(go.Scatter(
        x=data4.TIME,
        y=data4.ACC_Z,
        mode='lines',
        name=f'ACC Z where Y = -9.81',
    ), row=2, col=2)

    figure.add_trace(go.Scatter(
        x=data5.TIME,
        y=data5.ACC_X,
        mode='lines',
        name=f'ACC X where Z = 9.81',
    ), row=3, col=1)

    figure.add_trace(go.Scatter(
        x=data5.TIME,
        y=data5.ACC_Y,
        mode='lines',
        name=f'ACC Y where Z = 9.81',
    ), row=3, col=1)

    figure.add_trace(go.Scatter(
        x=data5.TIME,
        y=data5.ACC_Z,
        mode='lines',
        name=f'ACC Z where Z = 9.81',
    ), row=3, col=1)

    figure.add_trace(go.Scatter(
        x=data6.TIME,
        y=data6.ACC_X,
        mode='lines',
        name=f'ACC X where Z = -9.81',
    ), row=3, col=2)

    figure.add_trace(go.Scatter(
        x=data6.TIME,
        y=data6.ACC_Y,
        mode='lines',
        name=f'ACC Y where Z = -9.81',
    ), row=3, col=2)

    figure.add_trace(go.Scatter(
        x=data6.TIME,
        y=data6.ACC_Z,
        mode='lines',
        name=f'ACC Z where Z = -9.81',
    ), row=3, col=2)
    return figure


def loadFile(file, start=-1, ende=-1):
    Position_1 = pd.read_csv(file, sep=';')
    if start == -1:
        return Position_1[:ende]
    elif ende == -1:
        return Position_1[start:]
    elif ende == -1 and start == -1:
        return Position_1
    else:
        return Position_1[start:ende]


def accvelposplot(data):
    fig = make_subplots(3, 1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.ACC_X,
        mode='lines',
        name=f'ACC X',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.ACC_Y,
        mode='lines',
        name=f'ACC Y',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.ACC_Z,
        mode='lines',
        name=f'ACC Z',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.V_X,
        mode='lines',
        name=f'VX',
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.V_Y,
        mode='lines',
        name=f'VY',
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.V_Z,
        mode='lines',
        name=f'VZ',
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.P_X,
        mode='lines',
        name=f'PX',
    ), row=3, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.P_Y,
        mode='lines',
        name=f'PY',
    ), row=3, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.P_Z,
        mode='lines',
        name=f'PZ',
    ), row=3, col=1)

    return fig


def applyParams(data, params, offs):
    x = data['ACC_X']
    y = data['ACC_Y']
    z = data['ACC_Z']

    trueX = (x - offs[0]) / params[0][0]
    trueY = (y - offs[1]) / params[1][1]
    trueZ = (z - offs[2]) / params[2][2]

    df = data.copy()

    df['ACC_X'] = pd.Series(trueX, index=x.index)
    df['ACC_Y'] = pd.Series(trueY, index=y.index)
    df['ACC_Z'] = pd.Series(trueZ, index=z.index)

    return df


def plotAngles(data):
    fig = make_subplots(2, 1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.GYR_X,
        mode='lines',
        name=f'GYR X',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.GYR_Y,
        mode='lines',
        name=f'GYR Y',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.GYR_Z,
        mode='lines',
        name=f'GYR Z',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.RA_X,
        mode='lines',
        name=f'RA X',
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.RA_Y,
        mode='lines',
        name=f'RA Y',
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=data.TIME,
        y=data.RA_Z,
        mode='lines',
        name=f'RA Z',
    ), row=2, col=1)

    return fig


def load_baro(file, skiprows):
    df = pd.read_csv(file, skiprows=skiprows, delimiter="|")
    pr = df[df['sensorName'] == 'LPS22H Barometer Sensor'].reset_index(drop=True)
    del pr['statusId']
    del pr['sensorName']

    pr = pr[['timestamp', 'value']]
    vals = pr['value']
    vals = vals.str.replace('[', '')
    vals = vals.str.replace(']', '')
    x = [float(i) for i in vals]
    pr['pressure'] = x

    del pr['value']
    pr['unix_time'] = pd.to_datetime(pr['timestamp'], unit='ms')
    pr['unix_time'] -= timedelta(days=1, hours=21, minutes=36, seconds=41.235000)
    pr = pr.set_index(pr.unix_time)
    del pr['unix_time']
    pr = pr.resample('100ms').mean().interpolate()
    pr['unix_time'] = pr.index
    pr.reset_index(drop=True, inplace=True)

    del pr['timestamp']
    return pr


if __name__ == '__main__':
    gps = load_data('data/r10.pos', 23)
    print(gps)
    baro = load_baro("data/log.txt", 2)
    print(baro)

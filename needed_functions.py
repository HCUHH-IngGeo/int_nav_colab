import pandas as pd
import plotly.graph_objects as go
from pyproj import Proj, transform
from sympy import re, sqrt
from colour import Color

mapbox_access_token = "pk.eyJ1IjoianN0ODkiLCJhIjoiY2sybjN0b2w2MG1tYjNjcGV5eGtlYjg5NyJ9.A9sSZPYzTYJL1YErdmwoKA"


def load_data(file, skiprows, transform_to_utm=True):
    # read the csv file and skip rows correspondig to inpunt
    df = pd.read_csv(file, skiprows=skiprows)
    # renaming some colums to get easier access on columns
    df.rename(columns={df.columns[0]: "GPST", df.columns[1]: "lat", df.columns[2]: "lon", df.columns[3]: 'alt',
                       df.columns[6]: 'sdn', df.columns[7]: 'sde'},
              inplace=True)
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

    return df


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
    time = [time.strftime("%Y/%m/%d, %H:%M:%S") for time in df.GPST]
    figure.add_trace(go.Scattermapbox(
        lat=df.lat,
        lon=df.lon,
        mode='markers+lines',
        name=f'{name}',
        line={'color': color_of_trajectory},
        customdata=time,
        hovertemplate=
        "<b>%{customdata} </b><br><br>" +
        "longitude: %{lon}<br>" +
        "latitude: %{lat}<br>"
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

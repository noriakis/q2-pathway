
var path = require('path');
var webpack = require('webpack');

module.exports = {
  entry: {
    app: './src/main.js',
    vendor: ['d3']
  },
  plugins: [
    new webpack.optimize.CommonsChunkPlugin('vendor', 'dist/vendor.bundle.js'),
    new webpack.optimize.UglifyJsPlugin({
      compress: { warnings: false }
    }),
    new webpack.NoErrorsPlugin(),
  ],
  output: {
    path: __dirname,
    filename: 'dist/bundle.js'
  },
  module: {
    loaders: [
      {
        loader: 'babel-loader',
        exclude: /node_modules/,
        query: {
          presets: ['es2015']
        }
      }
    ]
  },
};

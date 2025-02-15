// pages/detail_edit/detail_edit.js
Page({

  /**
   * 页面的初始数据
   */
  data: {
    values:{

    }
  },

  /**
   * 生命周期函数--监听页面加载
   */
  onLoad:function(options) {
    if (options.item) {
      const receivedItem = JSON.parse(options.item);
      this.setData({
        values: receivedItem
      });
      console.log(this.data.values);
    }
  },
  updateValues: function(e) {
    const key = e.detail.key;
    const value = e.detail.value;
    this.setData({
      [`values.${key}`]: value
    });
    console.log(this.data.values);
  },
  /**
   * 生命周期函数--监听页面初次渲染完成
   */
  onReady() {

  },

  /**
   * 生命周期函数--监听页面显示
   */
  onShow() {

  },

  /**
   * 生命周期函数--监听页面隐藏
   */
  onHide() {

  },

  /**
   * 生命周期函数--监听页面卸载
   */
  onUnload() {

  },

  /**
   * 页面相关事件处理函数--监听用户下拉动作
   */
  onPullDownRefresh() {

  },

  /**
   * 页面上拉触底事件的处理函数
   */
  onReachBottom() {

  },

  /**
   * 用户点击右上角分享
   */
  onShareAppMessage() {

  }
})